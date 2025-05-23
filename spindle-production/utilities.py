"""Utilities module."""

import json
import os

import env_set
import networkx as nx
import numpy as np
import stk
import stko
from rdkit.Chem import AllChem as rdkit
from scipy.spatial.distance import euclidean
from stk import get_acute_vector, get_plane_normal


class MissingSettingError(Exception): ...


def read_lib(lib_file):
    """Read lib file."""
    with open(lib_file, "rb") as f:
        return json.load(f)


def get_atom_distance(molecule, atom1_id, atom2_id):
    """Return the distance between atom1 and atom2."""
    position_matrix = molecule.get_position_matrix()

    return float(
        euclidean(u=position_matrix[atom1_id], v=position_matrix[atom2_id])
    )


def reorient_linker(molecule):
    target_coords = (
        np.array([1, 1, 0]),
        np.array([1, -1, 0]),
        np.array([-1, -1, 0]),
        np.array([-1, 1, 0]),
    )
    centroid_pos = np.array([0, 0, 0])
    molecule = molecule.with_centroid(
        position=centroid_pos,
        atom_ids=molecule.get_placer_ids(),
    )

    edge_centroid = sum(target_coords) / len(target_coords)
    edge_normal = get_acute_vector(
        reference=edge_centroid,
        vector=get_plane_normal(
            points=np.array(target_coords),
        ),
    )

    fg_bonder_centroid = molecule.get_centroid(
        atom_ids=next(molecule.get_functional_groups()).get_placer_ids(),
    )
    edge_position = target_coords[0]
    molecule = molecule.with_rotation_to_minimize_angle(
        start=fg_bonder_centroid - centroid_pos,
        target=edge_position - edge_centroid,
        axis=edge_normal,
        origin=centroid_pos,
    )

    # Flatten wrt to xy plane.
    core_centroid = molecule.get_centroid(
        atom_ids=molecule.get_core_atom_ids(),
    )
    normal = molecule.get_plane_normal(
        atom_ids=molecule.get_placer_ids(),
    )
    normal = get_acute_vector(
        reference=core_centroid - centroid_pos,
        vector=normal,
    )
    molecule = molecule.with_rotation_between_vectors(
        start=normal,
        target=[0, 0, 1],
        origin=centroid_pos,
    )

    # Align long axis of molecule (defined by deleter atoms) with
    # y axis.
    long_axis_vector = molecule.get_long_axis()
    edge_centroid = sum(target_coords) / len(target_coords)
    edge_normal = get_acute_vector(
        reference=edge_centroid,
        vector=get_plane_normal(
            points=np.array(target_coords),
        ),
    )
    molecule = molecule.with_rotation_to_minimize_angle(
        start=long_axis_vector,
        target=[1, 0, 0],
        axis=edge_normal,
        origin=centroid_pos,
    )
    return molecule


def get_organic_linkers(cage, metal_atom_nos, output_dir, file_prefix=None):
    """Extract a list of organic linker .Molecules from a cage.

    Parameters
    ----------
    cage : :class:`stk.Molecule`
        Molecule to get the organic linkers from.

    metal_atom_nos : :class:`iterable` of :class:`int`
        The atomic number of metal atoms to remove from structure.

    file_prefix : :class:`str`, optional
        Prefix to file name of each output ligand structure.
        Eventual file name is:
        "file_prefix"{number of atoms}_{idx}_{i}.mol
        Where `idx` determines if a molecule is unique by smiles.

    Returns:
    -------
    org_lig : :class:`dict` of :class:`stk.BuildingBlock`
        Dictionary of building blocks where the key is the file name,
        and the value is the stk building block.

    smiles_keys : :class:`dict` of :class:`int`
        Key is the linker smiles, value is the idx of that smiles.

    """
    org_lig = {}

    # Produce a graph from the cage that does not include metals.
    cage_g = nx.Graph()
    atom_ids_in_G = set()
    for atom in cage.get_atoms():
        if atom.get_atomic_number() in metal_atom_nos:
            continue
        cage_g.add_node(atom)
        atom_ids_in_G.add(atom.get_id())

    # Add edges.
    for bond in cage.get_bonds():
        a1id = bond.get_atom1().get_id()
        a2id = bond.get_atom2().get_id()
        if a1id in atom_ids_in_G and a2id in atom_ids_in_G:
            cage_g.add_edge(bond.get_atom1(), bond.get_atom2())

    # Get disconnected subgraphs as molecules.
    # Sort and sort atom ids to ensure molecules are read by RDKIT
    # correctly.
    connected_graphs = [
        sorted(subgraph, key=lambda a: a.get_id())
        for subgraph in sorted(nx.connected_components(cage_g))
    ]
    smiles_keys = {}
    for i, cg in enumerate(connected_graphs):
        # Get atoms from nodes.
        atoms = list(cg)
        atom_ids = [i.get_id() for i in atoms]
        cage.write("temporary_linker.mol", atom_ids=atom_ids)
        temporary_linker = stk.BuildingBlock.init_from_file(
            "temporary_linker.mol"
        ).with_canonical_atom_ordering()
        smiles_key = stk.Smiles().get_key(temporary_linker)
        if smiles_key not in smiles_keys:
            smiles_keys[smiles_key] = len(smiles_keys.values()) + 1
        idx = smiles_keys[smiles_key]
        sgt = str(len(atoms))
        # Write to mol file.
        if file_prefix is None:
            filename_ = f"organic_linker_s{sgt}_{idx}_{i}.mol"
        else:
            filename_ = f"{file_prefix}{sgt}_{idx}_{i}.mol"

        org_lig[filename_] = temporary_linker
        os.system("rm temporary_linker.mol")
        # Rewrite to fix atom ids.
        org_lig[filename_].write(output_dir / filename_)
        org_lig[filename_] = stk.BuildingBlock.init_from_file(
            str(output_dir / filename_)
        )

    return org_lig, smiles_keys


def get_lowest_energy_conformers(
    org_ligs,
    smiles_keys,
    file_prefix,
    settings,
    output_dir,
):
    """Determine the lowest energy conformer of cage organic linkers.

    Will do multiple if there are multiple types.

    Uses previously run crest outputs (neew crest does not work with old xtb!).

    Parameters
    ----------
    org_ligs : :class:`dict` of :class:`stk.BuildingBlock`
        Dictionary of building blocks where the key is the file name,
        and the value is the stk building block.

    smiles_keys : :class:`dict` of :class:`int`
        Key is the linker smiles, value is the idx of that smiles.

    file_prefix : :class:`str`
        Prefix to file name of each output ligand structure.
        Eventual file name is:
        "file_prefix"{number of atoms}_{idx}_{i}.mol
        Where `idx` determines if a molecule is unique by smiles.

    """
    for lig in org_ligs:
        stk_lig = org_ligs[lig]
        smiles_key = stk.Smiles().get_key(stk_lig)
        idx = smiles_keys[smiles_key]
        sgt = str(stk_lig.get_num_atoms())
        final_filename_ = output_dir / f"{file_prefix}{sgt}_{idx}_opt.mol"
        ligand_name_ = file_prefix.split("_sg")[0].split("_")[2:]
        if len(ligand_name_) > 1:
            ligand_name_ = "_".join(ligand_name_)
        else:
            ligand_name_ = ligand_name_[0]

        if not final_filename_.exists():
            msg = f"File {final_filename_} does not exist."
            raise RuntimeError(msg)


def get_dihedral(pt1, pt2, pt3, pt4):
    """Calculate the dihedral between four points.

    Uses Praxeolitic formula --> 1 sqrt, 1 cross product

    Output in range (-pi to pi).

    From: https://stackoverflow.com/questions/20305272/
    dihedral-torsion-angle-from-four-points-in-cartesian-
    coordinates-in-python
    (new_dihedral(p))

    """
    p0 = np.asarray(pt1)
    p1 = np.asarray(pt2)
    p2 = np.asarray(pt3)
    p3 = np.asarray(pt4)

    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def get_query_atom_ids(query, rdkit_mol):
    """Yield the ids of atoms in `rdkit_mol` which match `query`.

    Multiple substructures in `rdkit_mol` can match `query` and
    therefore each set is yielded as a group.

    Parameters
    ----------
    query : :class:`str`
        A SMARTS string used to query atoms.

    rdkit_mol : :class:`rdkit.Mol`
        A molecule whose atoms should be queried.

    Yields:
    ------
    :class:`tuple` of :class:`int`
        The ids of atoms in `molecule` which match `query`.

    """
    rdkit.SanitizeMol(rdkit_mol)
    yield from rdkit_mol.GetSubstructMatches(
        query=rdkit.MolFromSmarts(query),
    )


def calculate_abs_imine_torsions(org_ligs, smarts=None):
    """Calculate the imine torsion of all ligands in the cage."""
    if smarts is None:
        # C-N=C(H)-C(X)-X, where X != H.
        smarts = "[#6]-[#7X2]=[#6X3H1]-[#6X3!H1]"

    torsions = {}
    # Iterate over ligands.
    for lig in org_ligs:
        stk_lig = org_ligs[lig]
        # Find torsions.
        rdkit_mol = stk_lig.to_rdkit_mol()
        query_ids = get_query_atom_ids(smarts, rdkit_mol)
        # Calculate torsional angle for all imines.
        torsion_list = []
        for atom_ids in query_ids:
            torsion = get_dihedral(
                pt1=next(iter(stk_lig.get_atomic_positions(atom_ids[0]))),
                pt2=next(iter(stk_lig.get_atomic_positions(atom_ids[1]))),
                pt3=next(iter(stk_lig.get_atomic_positions(atom_ids[2]))),
                pt4=next(iter(stk_lig.get_atomic_positions(atom_ids[3]))),
            )
            torsion_list.append(abs(torsion))

        # Degrees
        torsions[lig] = torsion_list

    return torsions


def read_gfnx2xtb_eyfile(file):
    """Read the energy (kJ/mol from GFN2-xTB) from a .ey file."""
    with open(file) as f:
        lines = f.readlines()
        ey = float(lines[0].rstrip())

    return ey * 2625.5


def calculate_energy(
    name,
    mol,
    ey_file,
    output_dir,
    xtb_path=None,
    charge=0,
    no_unpaired_e=0,
    solvent=None,
):
    """Calculate GFN-xTB energy of molecule."""
    if xtb_path is None:
        xtb_path = env_set.xtb_path()

    print(f"....getting energy of {name}")
    xtb_energy = stko.XTBEnergy(
        xtb_path=env_set.xtb_path(),
        output_dir=output_dir / f"{name}_ey",
        num_cores=6,
        charge=charge,
        num_unpaired_electrons=no_unpaired_e,
        electronic_temperature=300,
        unlimited_memory=True,
        calculate_free_energy=False,
        solvent=solvent,
    )
    energy = xtb_energy.get_energy(mol)

    with open(ey_file, "w") as f:
        f.write(str(energy))


def calculate_ligand_SE(
    org_ligs,
    smiles_keys,
    output_json,
    output_dir,
    file_prefix=None,
    solvent=None,
):
    """Calculate the strain energy of each ligand in the cage.

    Parameters
    ----------
    org_lig : :class:`dict` of :class:`stk.BuildingBlock`
        Dictionary of building blocks where the key is the file name,
        and the value is the stk building block.

    smiles_keys : :class:`dict` of :class:`int`
        Key is the linker smiles, value is the idx of that smiles.

    output_json : :class:`str`
        File name to save output to to avoid reruns.

    file_prefix : :class:`str`, optional
        Prefix to file name of each output ligand structure.
        Eventual file name is:
        "file_prefix"{number of atoms}_{idx}_{i}.mol
        Where `idx` determines if a molecule is unique by smiles.

    solvent: :class:`str`
        None if gas phase, otherwise a string matching a solvent
        model available in xtb.

    Returns:
    -------
    strain_energies : :class:`dict`
        Strain energies for each ligand.

    """
    # Check if output file exists.
    if not os.path.exists(output_json):
        strain_energies = {}
        # Iterate over ligands.
        for lig in org_ligs:
            stk_lig = org_ligs[lig]
            ey_file = output_dir / lig.replace("mol", "ey")
            smiles_key = stk.Smiles().get_key(stk_lig)
            idx = smiles_keys[smiles_key]
            sgt = str(stk_lig.get_num_atoms())
            # Get optimized ligand name that excludes any cage
            # information.
            if file_prefix is None:
                filename_ = output_dir / f"organic_linker_s{sgt}_{idx}_opt.mol"
                opt_lig_ey = output_dir / f"organic_linker_s{sgt}_{idx}_opt.ey"
                opt_lig_n = f"organic_linker_s{sgt}_{idx}_opt"
            else:
                filename_ = output_dir / f"{file_prefix}{sgt}_{idx}_opt.mol"
                opt_lig_ey = output_dir / f"{file_prefix}{sgt}_{idx}_opt.ey"
                opt_lig_n = f"{file_prefix}{sgt}_{idx}_opt"

            # Calculate energy of extracted ligand.
            if not os.path.exists(ey_file):
                calculate_energy(
                    name=lig.replace(".mol", ""),
                    mol=stk_lig,
                    ey_file=ey_file,
                    output_dir=output_dir,
                    solvent=solvent,
                )
            # Read energy.
            # kJ/mol.
            E_extracted = read_gfnx2xtb_eyfile(ey_file)

            # Calculate energy of optimised ligand.
            # Load in lowest energy conformer.
            opt_mol = stk.BuildingBlock.init_from_file(str(filename_))
            if not os.path.exists(opt_lig_ey):
                calculate_energy(
                    name=opt_lig_n,
                    mol=opt_mol,
                    ey_file=opt_lig_ey,
                    solvent=solvent,
                    output_dir=output_dir,
                )
            # Read energy.
            # kJ/mol.
            print(opt_lig_ey, filename_, opt_lig_n)
            E_free = read_gfnx2xtb_eyfile(opt_lig_ey)
            # Add to list the strain energy:
            # (E(extracted) - E(optimised/free))
            lse = E_extracted - E_free
            # kJ/mol.
            strain_energies[lig] = lse

        # Write data.
        with open(output_json, "w") as f:
            json.dump(strain_energies, f)

    # Get data.
    with open(output_json) as f:
        strain_energies = json.load(f)

    return strain_energies


def calculate_metal_ligand_distance(
    mol,
    metal_atomic_number,
    ligand_atomic_number,
):
    """Calculate all bond lengths in mol between metal and ligand atoms.

    Parameters
    ----------
    mol : :class:`stk.ConstructedMolecule`
        stk molecule to analyse.

    metal_atomic_number : :class:`int`
        Element number of metal atom.

    ligand_atomic_number : :class:`int`
        Element number of atoms bonded to metal.

    Returns:
    -------
    bond_lengths : :class:`list`
        Bond lengths (in Angstrom).

    """
    bond_lengths = []
    # Calculate bond lengths.
    for bond in mol.get_bonds():
        atom1_id = bond.get_atom1().get_id()
        atom2_id = bond.get_atom2().get_id()
        atom1_an = bond.get_atom1().get_atomic_number()
        atom2_an = bond.get_atom2().get_atomic_number()
        chk1 = (
            atom1_an == metal_atomic_number
            and atom2_an == ligand_atomic_number
        )
        chk2 = (
            atom1_an == ligand_atomic_number
            and atom2_an == metal_atomic_number
        )
        if chk1 or chk2:
            bond_lengths.append(
                get_atom_distance(
                    molecule=mol,
                    atom1_id=atom1_id,
                    atom2_id=atom2_id,
                )
            )
    return bond_lengths
