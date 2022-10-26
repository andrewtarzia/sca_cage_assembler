#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Utilities module.

Author: Andrew Tarzia

Date Created: 15 Mar 2020
"""

import subprocess as sp
import numpy as np
import networkx as nx
import os
import json
from itertools import combinations, permutations
from rdkit.Chem import AllChem as rdkit
from rdkit.Chem import Draw
import pymatgen.core as pmg
from pymatgen.analysis.local_env import (
    LocalStructOrderParams,
)
from scipy.spatial.distance import euclidean
import stk
from stk import (
    get_acute_vector,
    get_plane_normal,
)
import stko

import env_set


class MissingSettingError(Exception):
    ...


def get_lowest_energy_conformer(name, mol, conf_dir, settings):
    """
    Get lowest energy conformer of molecule.

    Method:
        1) squick CREST conformer search
        2) xTB `opt_level` optimisation of lowest energy conformer
        3) save file

    """

    # Check for missing settings.
    req_settings = [
        'final_opt_level', 'conf_opt_level', 'charge', 'no_unpaired_e',
        'max_runs', 'calc_hessian', 'solvent', 'nc', 'crest_exec',
        'etemp', 'keepdir', 'cross', 'md_len', 'ewin', 'speed_setting'
    ]
    for i in req_settings:
        if i not in settings:
            raise MissingSettingError(
                f'Settings missing {i}. Has {settings.keys()}.'
            )

    low_e_conf = crest_conformer_search(
        molecule=mol,
        output_dir=conf_dir,
        gfn_version=2,
        nc=settings['nc'],
        opt_level=settings['conf_opt_level'],
        charge=settings['charge'],
        etemp=settings['etemp'],
        no_unpaired_e=settings['no_unpaired_e'],
        keepdir=settings['keepdir'],
        cross=settings['cross'],
        md_len=settings['md_len'],
        ewin=settings['ewin'],
        speed_setting=settings['speed_setting'],
        solvent=settings['solvent'],
    )

    # Save lowest energy conformer.
    low_e_conf.write(os.path.join(conf_dir, 'low_e_unopt.mol'))

    # Optimize lowest energy conformer at opt_level.
    low_e_conf = optimize_conformer(
        name=f'{name}_low_e_opt',
        mol=low_e_conf,
        opt_level=settings['final_opt_level'],
        charge=settings['charge'],
        no_unpaired_e=settings['no_unpaired_e'],
        max_runs=settings['max_runs'],
        calc_hessian=settings['calc_hessian'],
        solvent=settings['solvent']
    )
    low_e_conf.write(os.path.join(conf_dir, 'low_e_opt.mol'))

    # Return molecule.
    return low_e_conf


def get_lowest_energy_conformers(
    org_ligs,
    smiles_keys,
    file_prefix,
    settings,
):
    """
    Determine the lowest energy conformer of cage organic linkers.

    Will do multiple if there are multiple types.

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
        final_filename_ = f'{file_prefix}{sgt}_{idx}_opt.mol'
        ligand_name_ = file_prefix.split('_sg')[0].split('_')[2:]
        if len(ligand_name_) > 1:
            ligand_name_ = '_'.join(ligand_name_)
        else:
            ligand_name_ = ligand_name_[0]

        if os.path.exists(final_filename_):
            continue
        else:
            print(
                '......calculating lowest energy conformer of '
                f'{ligand_name_}'
            )
            # Get low energy conformer using CREST.
            conf_dir = f'{ligand_name_}_xtbcrest_confs'
            if not os.path.exists(conf_dir):
                os.mkdir(conf_dir)
            low_e_conf = get_lowest_energy_conformer(
                name=ligand_name_,
                mol=stk_lig,
                conf_dir=conf_dir,
                settings=env_set.crest_conformer_settings(
                    solvent=None
                ),
            )
            # In cage specific case, we want this optimised with solvent.
            low_e_conf = optimize_conformer(
                name=ligand_name_+'low_e_opt',
                mol=low_e_conf,
                opt_level=settings['final_opt_level'],
                charge=settings['charge'],
                no_unpaired_e=settings['no_unpaired_e'],
                max_runs=settings['max_runs'],
                calc_hessian=settings['calc_hessian'],
                solvent=settings['solvent']
            )
            low_e_conf.write(final_filename_)


def optimize_conformer(
    name,
    mol,
    opt_level='extreme',
    charge=0,
    no_unpaired_e=0,
    max_runs=1,
    calc_hessian=False,
    solvent=None
):
    """
    Run simple GFN-xTB optimisation of molecule.

    """

    print(f'....optimizing {name}')
    xtb_opt = stko.XTB(
        xtb_path=env_set.xtb_path(),
        output_dir=f'{name}_opt',
        gfn_version=2,
        num_cores=6,
        opt_level=opt_level,
        charge=charge,
        num_unpaired_electrons=no_unpaired_e,
        max_runs=max_runs,
        calculate_hessian=calc_hessian,
        unlimited_memory=True,
        solvent=solvent,
    )

    return xtb_opt.optimize(mol=mol)


def crest_conformer_search(
    molecule,
    output_dir,
    gfn_version,
    nc,
    opt_level,
    charge,
    keepdir,
    cross,
    etemp,
    no_unpaired_e,
    md_len=None,
    ewin=5,
    speed_setting=None,
    solvent=None,
):
    """
    Perform conformer search with GFN-2 and CREST.

    Parameters
    ----------

    Returns
    -------

    """

    print('..........doing GFN-2 CREST optimisation')
    xtb_crest = stko.XTBCREST(
        crest_path=env_set.crest_path(),
        xtb_path=env_set.xtb_path(),
        gfn_version=gfn_version,
        output_dir=output_dir,
        num_cores=nc,
        ewin=ewin,
        opt_level=opt_level,
        charge=charge,
        electronic_temperature=etemp,
        num_unpaired_electrons=no_unpaired_e,
        solvent=solvent,
        keepdir=keepdir,
        cross=cross,
        speed_setting=speed_setting,
        md_len=md_len,
        unlimited_memory=True,
    )
    molecule = xtb_crest.optimize(mol=molecule)

    return molecule


def convert_stk_to_pymatgen(stk_mol):
    """
    Convert stk.Molecule to pymatgen.Molecule.

    Parameters
    ----------
    stk_mol : :class:`stk.Molecule`
        Stk molecule to convert.

    Returns
    -------
    pmg_mol : :class:`pymatgen.Molecule`
        Corresponding pymatgen Molecule.

    """
    stk_mol.write('temp.xyz')
    pmg_mol = pmg.Molecule.from_file('temp.xyz')
    os.system('rm temp.xyz')

    return pmg_mol


def get_element_sites(molecule, atomic_no):
    """
    Get the index of sites in the molecule with the desired element.

    Parameters
    ----------
    molecule : :class:`pymatgen.Molecule`
        Molecule to get sites of.

    atomic_no : :class:`int`
        Atomic number of desired element.

    Returns
    -------
    idxs : :class:`list` of :class:`int`
        List of site indices.

    """

    idxs = []

    for i, atom in enumerate(molecule.sites):
        Z = atom.specie.Z
        if Z == atomic_no:
            idxs.append(i)

    return idxs


def calculate_sites_order_values(
    molecule,
    site_idxs,
    target_species_type=None,
    neigh_idxs=None
):
    """
    Calculate order parameters around metal centres.

    Parameters
    ----------
    molecule : :class:`pmg.Molecule` or :class:`pmg.Structure`
        Pymatgen (pmg) molecule/structure to analyse.

    site_idxs : :class:`list` of :class:`int`
        Atom ids of sites to calculate OP of.

    target_species_type : :class:`str`
        Target neighbour element to use in OP calculation.
        Defaults to :class:`NoneType` if no target species is known.

    neigh_idxs : :class:`list` of :class:`list` of :class:`int`
        Neighbours of each atom in site_idx. Ordering is important.
        Defaults to :class:`NoneType` for when using
        :class:`pmg.Structure` - i.e. a structure with a lattice.

    Returns
    -------
    results : :class:`dict`
        Dictionary of format
        site_idx: dict of order parameters
        {
            `oct`: :class:`float`,
            `sq_plan`: :class:`float`,
            `q2`: :class:`float`,
            `q4`: :class:`float`,
            `q6`: :class:`float`
        }.

    """

    results = {}

    if target_species_type is None:
        targ_species = None
    else:
        targ_species = pmg.Species(target_species_type)

    # Define local order parameters class based on desired types.
    types = [
        'oct',  # Octahedra OP.
        'sq_plan',  # Square planar envs.
        'q2',  # l=2 Steinhardt OP.
        'q4',  # l=4 Steinhardt OP.
        'q6',  # l=6 Steinhardt OP.
    ]
    loc_ops = LocalStructOrderParams(
        types=types,
    )
    if neigh_idxs is None:
        for site in site_idxs:
            site_results = loc_ops.get_order_parameters(
                structure=molecule,
                n=site,
                target_spec=[targ_species]
            )
            results[site] = {i: j for i, j in zip(types, site_results)}
    else:
        for site, neigh in zip(site_idxs, neigh_idxs):
            site_results = loc_ops.get_order_parameters(
                structure=molecule,
                n=site,
                indices_neighs=neigh,
                target_spec=targ_species
            )
            results[site] = {i: j for i, j in zip(types, site_results)}

    return results


def get_order_values(mol, metal, per_site=False):
    """
    Calculate order parameters around metal centres.

    Parameters
    ----------
    mol : :class:`stk.ConstructedMolecule`
        stk molecule to analyse.

    metal : :class:`int`
        Element number of metal atom.

    per_site : :class:`bool`
        Defaults to False. True if the OPs for each site are desired.

    Returns
    -------
    results : :class:`dict`
        Dictionary of order parameter max/mins/averages if `per_site`
        is False.

    """

    pmg_mol = convert_stk_to_pymatgen(stk_mol=mol)
    # Get sites of interest and their neighbours.
    sites = []
    neighs = []
    for atom in mol.get_atoms():
        if atom.get_atomic_number() == metal:
            sites.append(atom.get_id())
            bonds = [
                i
                for i in mol.get_bonds()
                if i.get_atom1().get_id() == atom.get_id()
                or i.get_atom2().get_id() == atom.get_id()
            ]
            a_neigh = []
            for b in bonds:
                if b.get_atom1().get_id() == atom.get_id():
                    a_neigh.append(b.get_atom2().get_id())
                elif b.get_atom2().get_id() == atom.get_id():
                    a_neigh.append(b.get_atom1().get_id())
            neighs.append(a_neigh)

    order_values = calculate_sites_order_values(
        molecule=pmg_mol,
        site_idxs=sites,
        neigh_idxs=neighs
    )

    if per_site:
        results = order_values
        return results
    else:
        # Get max, mins and averages of all OPs for the whole molecule.
        OPs = [order_values[i].keys() for i in order_values][0]
        OP_lists = {}
        for OP in OPs:
            OP_lists[OP] = [order_values[i][OP] for i in order_values]

        results = {
            # OP: (min, max, avg)
            i: {
                'min': min(OP_lists[i]),
                'max': max(OP_lists[i]),
                'avg': np.average(OP_lists[i])
            }
            for i in OP_lists
        }

        return results


def get_organic_linkers(cage, metal_atom_nos, file_prefix=None):
    """
    Extract a list of organic linker .Molecules from a cage.

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

    Returns
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
        cage.write(
            'temporary_linker.mol',
            atom_ids=atom_ids
        )
        temporary_linker = stk.BuildingBlock.init_from_file(
            'temporary_linker.mol'
        ).with_canonical_atom_ordering()
        smiles_key = stk.Smiles().get_key(temporary_linker)
        if smiles_key not in smiles_keys:
            smiles_keys[smiles_key] = len(smiles_keys.values())+1
        idx = smiles_keys[smiles_key]
        sgt = str(len(atoms))
        # Write to mol file.
        if file_prefix is None:
            filename_ = f'organic_linker_s{sgt}_{idx}_{i}.mol'
        else:
            filename_ = f'{file_prefix}{sgt}_{idx}_{i}.mol'

        org_lig[filename_] = temporary_linker
        os.system('rm temporary_linker.mol')
        # Rewrite to fix atom ids.
        org_lig[filename_].write(filename_)
        org_lig[filename_] = stk.BuildingBlock.init_from_file(
            filename_
        )

    return org_lig, smiles_keys


def draw_and_save_grid(
    mol_list,
    names,
    subImgSize,
    mol_per_row,
    filename
):
    """
    Draw RDKit molecules and save SVG.

    """
    img = Draw.MolsToGridImage(
        mol_list,
        molsPerRow=mol_per_row,
        subImgSize=subImgSize,
        legends=names,
        useSVG=True
    )
    save_svg(
        filename=filename,
        string=img
    )


def mol_list2grid(
    molecules,
    filename,
    mol_per_row,
    maxrows,
    subImgSize=(200, 200),
    names=None
):
    """
    Produce a grid of molecules in mol_list.

    molecules (list) - list of molecule SMILEs

    """

    if len(molecules) > mol_per_row * maxrows:
        # have to make multiple images
        new_mol_list = []
        new_names = []
        count = 1
        for i, mol in enumerate(molecules):
            new_mol_list.append(mol)
            if names is None:
                new_names = None
            else:
                new_names.append(names[i])
            # make image
            chk1 = len(new_mol_list) == mol_per_row * maxrows
            chk2 = i == len(molecules)-1
            if chk1 or chk2:
                draw_and_save_grid(
                    mol_list=new_mol_list,
                    mol_per_row=mol_per_row,
                    subImgSize=subImgSize,
                    names=new_names,
                    filename=f'{filename}_{count}.svg'
                )
                # img.save(filename + '_' + str(count) + '.png')
                new_mol_list = []
                new_names = []
                count += 1
    else:
        draw_and_save_grid(
            mol_list=molecules,
            mol_per_row=mol_per_row,
            subImgSize=subImgSize,
            names=names,
            filename=f'{filename}.svg'
        )


def save_svg(filename, string):
    """
    Save svg text to a file.

    """

    with open(filename, 'w') as f:
        f.write(string)


def read_gfnx2xtb_eyfile(file):
    """
    Read the energy (kJ/mol from GFN2-xTB) from a .ey file.

    """

    with open(file, 'r') as f:
        lines = f.readlines()
        ey = float(lines[0].rstrip())

    return ey*2625.5


def calculate_energy(
    name,
    mol,
    ey_file,
    xtb_path=None,
    charge=0,
    no_unpaired_e=0,
    solvent=None
):
    """
    Calculate GFN-xTB energy of molecule.

    """

    if xtb_path is None:
        xtb_path = env_set.xtb_path()

    print(f'....getting energy of {name}')
    xtb_energy = stko.XTBEnergy(
        xtb_path=env_set.xtb_path(),
        output_dir=f'{name}_ey',
        num_cores=6,
        charge=charge,
        num_unpaired_electrons=no_unpaired_e,
        electronic_temperature=300,
        unlimited_memory=True,
        calculate_free_energy=False,
        solvent=solvent,
    )
    energy = xtb_energy.get_energy(mol)

    with open(ey_file, 'w') as f:
        f.write(str(energy))


def calculate_ligand_SE(
    org_ligs,
    smiles_keys,
    output_json,
    file_prefix=None,
    solvent=None,
):
    """
    Calculate the strain energy of each ligand in the cage.

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

    Returns
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
            ey_file = lig.replace('mol', 'ey')
            smiles_key = stk.Smiles().get_key(stk_lig)
            idx = smiles_keys[smiles_key]
            sgt = str(stk_lig.get_num_atoms())
            # Get optimized ligand name that excludes any cage
            # information.
            if file_prefix is None:
                filename_ = f'organic_linker_s{sgt}_{idx}_opt.mol'
                opt_lig_ey = f'organic_linker_s{sgt}_{idx}_opt.ey'
                opt_lig_n = f'organic_linker_s{sgt}_{idx}_opt'
            else:
                filename_ = f'{file_prefix}{sgt}_{idx}_opt.mol'
                opt_lig_ey = f'{file_prefix}{sgt}_{idx}_opt.ey'
                opt_lig_n = f'{file_prefix}{sgt}_{idx}_opt'

            # Calculate energy of extracted ligand.
            if not os.path.exists(ey_file):
                calculate_energy(
                    name=lig.replace('.mol', ''),
                    mol=stk_lig,
                    ey_file=ey_file,
                    solvent=solvent,
                )
            # Read energy.
            # kJ/mol.
            E_extracted = read_gfnx2xtb_eyfile(ey_file)

            # Calculate energy of optimised ligand.
            # Load in lowest energy conformer.
            opt_mol = stk.BuildingBlock.init_from_file(
                filename_
            )
            if not os.path.exists(opt_lig_ey):
                calculate_energy(
                    name=opt_lig_n,
                    mol=opt_mol,
                    ey_file=opt_lig_ey,
                    solvent=solvent,
                )
            # Read energy.
            # kJ/mol.
            print(opt_lig_ey, filename_, opt_lig_n)
            input()
            E_free = read_gfnx2xtb_eyfile(opt_lig_ey)
            # Add to list the strain energy:
            # (E(extracted) - E(optimised/free))
            lse = E_extracted - E_free
            # kJ/mol.
            strain_energies[lig] = lse

        # Write data.
        with open(output_json, 'w') as f:
            json.dump(strain_energies, f)

    # Get data.
    with open(output_json, 'r') as f:
        strain_energies = json.load(f)

    return strain_energies


def get_query_atom_ids(query, rdkit_mol):
    """
    Yield the ids of atoms in `rdkit_mol` which match `query`.

    Multiple substructures in `rdkit_mol` can match `query` and
    therefore each set is yielded as a group.

    Parameters
    ----------
    query : :class:`str`
        A SMARTS string used to query atoms.

    rdkit_mol : :class:`rdkit.Mol`
        A molecule whose atoms should be queried.

    Yields
    ------
    :class:`tuple` of :class:`int`
        The ids of atoms in `molecule` which match `query`.

    """

    rdkit.SanitizeMol(rdkit_mol)
    yield from rdkit_mol.GetSubstructMatches(
        query=rdkit.MolFromSmarts(query),
    )


def get_dihedral(pt1, pt2, pt3, pt4):
    """
    Calculate the dihedral between four points.

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


def calculate_abs_imine_torsions(org_ligs, smarts=None):
    """
    Calculate the imine torsion of all ligands in the cage.

    """

    if smarts is None:
        # C-N=C(H)-C(X)-X, where X != H.
        smarts = '[#6]-[#7X2]=[#6X3H1]-[#6X3!H1]'

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
                pt1=tuple(
                    stk_lig.get_atomic_positions(atom_ids[0])
                )[0],
                pt2=tuple(
                    stk_lig.get_atomic_positions(atom_ids[1])
                )[0],
                pt3=tuple(
                    stk_lig.get_atomic_positions(atom_ids[2])
                )[0],
                pt4=tuple(
                    stk_lig.get_atomic_positions(atom_ids[3])
                )[0]
            )
            torsion_list.append(abs(torsion))

        # Degrees
        torsions[lig] = torsion_list

    return torsions


def calculate_ligand_planarities(org_ligs):
    """
    Calculate the planarity of all ligands.

    """

    planarities = {}
    # Iterate over ligands.
    for lig in org_ligs:
        stk_lig = org_ligs[lig]

        # Calculate planarity.
        planarity = calculate_molecule_planarity(stk_lig)

        # Angstroms.
        planarities[lig] = planarity

    return planarities


def get_atom_distance(molecule, atom1_id, atom2_id):
    """
    Return the distance between atom1 and atom2.

    Parameters
    ----------
    molecule : :class:`stk.Molecule`

    atom1_id : :class:`int`
        The id of atom1.

    atom2_id : :class:`int`
        The id of atom2.

    Returns
    -------
    :class:`float`
        The euclidean distance between two atoms.

    """

    position_matrix = molecule.get_position_matrix()

    distance = euclidean(
        u=position_matrix[atom1_id],
        v=position_matrix[atom2_id]
    )

    return float(distance)


def shortest_distance_to_plane(plane, point):
    """
    Calculate the perpendicular distance beween a point and a plane.

    """

    top = abs(
        plane[0]*point[0] + plane[1]*point[1] +
        plane[2]*point[2] - plane[3]
    )
    bottom = np.sqrt(plane[0]**2 + plane[1]**2 + plane[2]**2)
    distance = top / bottom
    return distance


def calculate_molecule_planarity(mol, plane_ids=None, atom_ids=None):
    """
    Calculate planrity of molecule as sum of deviation from plane.

    Returns sum in Angstrom.

    Arguments
    ---------
    mol : :class:`stk.Molecule`
        Molecule.

    plane_ids : iterable of :class:`int`, optional
        Atom ids to use to define plane. Defaults to all atom ids.

    atom_ids : iterable of :class:`int`, optional
        Atom ids to calculate deviation for. Defaults to all atom ids.

    """

    if plane_ids is None:
        plane_ids = list(range(len(list(mol.get_atoms()))))
    else:
        plane_ids = list(plane_ids)

    if atom_ids is None:
        atom_ids = list(range(len(list(mol.get_atoms()))))
    else:
        atom_ids = list(atom_ids)

    centroid = mol.get_centroid(atom_ids=plane_ids)
    normal = mol.get_plane_normal(atom_ids=plane_ids)
    # Plane of equation ax + by + cz = d.
    atom_plane = np.append(normal, np.sum(normal*centroid))
    # Define the plane deviation as the sum of the distance of all
    # atoms from the plane defined by all atoms.
    plane_dev = sum([
        shortest_distance_to_plane(
            atom_plane,
            tuple(mol.get_atomic_positions(atom_ids=i.get_id()), )[0]
        )
        for i in mol.get_atoms() if i.get_id() in atom_ids
    ])

    return plane_dev


def unit_vector(vector):
    """
    Returns the unit vector of the vector.

    https://stackoverflow.com/questions/2827393/
    angles-between-two-n-dimensional-vectors-in-python/
    13849249#13849249
    """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2, normal=None):
    """
    Returns the angle in radians between vectors 'v1' and 'v2'::

        >>> angle_between((1, 0, 0), (0, 1, 0))
        1.5707963267948966
        >>> angle_between((1, 0, 0), (1, 0, 0))
        0.0
        >>> angle_between((1, 0, 0), (-1, 0, 0))
        3.141592653589793

    https://stackoverflow.com/questions/2827393/
    angles-between-two-n-dimensional-vectors-in-python/
    13849249#13849249

    If normal is given, the angle polarity is determined using the
    cross product of the two vectors.

    """

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    if normal is not None:
        # Get normal vector and cross product to determine sign.
        cross = np.cross(v1_u, v2_u)
        if np.dot(normal, cross) < 0:
            angle = -angle
    return angle


def read_lib(lib_file):
    """
    Read lib file.

    Returns dictionary.

    """

    print(f'reading {lib_file}')
    with open(lib_file, 'rb') as f:
        lib = json.load(f)

    return lib


def calculate_binding_AR(mol, atom_ids=None):
    """
    Calculate ligand aspect ratio based on binder positions.

    Defined as:
        Average ratio of the two shortest binding atom-binding atom
        distances eminating from each binding atom in the molecule.

    Using the atom_ids argument, you can define a different set of four
    atoms to define the anisotropy.

    """

    if atom_ids is None:
        if mol.get_num_functional_groups() != 4:
            return None
        target_atom_ids = [
            list(fg.get_bonder_ids())
            for fg in mol.get_functional_groups()
        ]
    else:
        if len(atom_ids) != 4:
            raise ValueError('Requires exactly four target atom ids.')
        target_atom_ids = atom_ids

    atom_dists = sorted(
        [
            (idx1, idx2, get_atom_distance(mol, idx1, idx2))
            for idx1, idx2 in combinations(target_atom_ids, r=2)
        ],
        key=lambda a: a[2]
    )

    far_binder_pair = (
        atom_dists[-1][0],
        atom_dists[-1][1]
    )
    ARs = []
    for fg_id in far_binder_pair:
        ds = sorted([
            i[2] for i in atom_dists
            if fg_id in (i[0], i[1])
        ])
        AR = ds[1]/min(ds)
        ARs.append(AR)

    ligand_AR = sum(ARs)/len(ARs)
    return ligand_AR


def calculate_binding_ABS(mol, atom_ids=None):
    """
    Calculate ligand binder differences based on binder positions.

    Defined as:
        Average diffference of the two shortest binding atom-binding
        atom distances eminating from each binding atom in the
        molecule.

    Using the atom_ids argument, you can define a different set of four
    atoms to define the anisotropy.

    """

    if atom_ids is None:
        if mol.get_num_functional_groups() != 4:
            return None
        target_atom_ids = [
            list(fg.get_bonder_ids())
            for fg in mol.get_functional_groups()
        ]
    else:
        if len(atom_ids) != 4:
            raise ValueError('Requires exactly four target atom ids.')
        target_atom_ids = atom_ids

    atom_dists = sorted(
        [
            (idx1, idx2, get_atom_distance(mol, idx1, idx2))
            for idx1, idx2 in combinations(target_atom_ids, r=2)
        ],
        key=lambda a: a[2]
    )

    far_binder_pair = (
        atom_dists[-1][0],
        atom_dists[-1][1]
    )
    ABSs = []
    for fg_id in far_binder_pair:
        ds = sorted([
            i[2] for i in atom_dists
            if fg_id in (i[0], i[1])
        ])
        AB = ds[1]-ds[0]
        ABSs.append(AB)

    ligand_AB = sum(ABSs)/len(ABSs)
    return ligand_AB


def calculate_metal_ligand_distance(
    mol,
    metal_atomic_number,
    ligand_atomic_number
):
    """
    Calculate all bond lengths in mol between metal and ligand atoms.

    Parameters
    ----------
    mol : :class:`stk.ConstructedMolecule`
        stk molecule to analyse.

    metal_atomic_number : :class:`int`
        Element number of metal atom.

    ligand_atomic_number : :class:`int`
        Element number of atoms bonded to metal.

    Returns
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


def calculate_ideal_pore_size(mol):
    """
    Calculate ideal pore size based on ligand deleter positions.

    Defined as:
        Avaerage distance between deleter atoms along short axis.

    """

    if mol.get_num_functional_groups() != 4:
        return None

    deleter_atom_ids = [
        list(fg.get_deleter_ids())
        for fg in mol.get_functional_groups()
    ]
    deleter_atom_dists = sorted(
        [
            (idx1, idx2, get_atom_distance(mol, idx1, idx2))
            for idx1, idx2 in combinations(deleter_atom_ids, r=2)
        ],
        key=lambda a: a[2]
    )

    far_binder_pair = (
        deleter_atom_dists[-1][0],
        deleter_atom_dists[-1][1]
    )
    short_distances = []
    for fg_id in far_binder_pair:
        ds = sorted([
            i[2] for i in deleter_atom_dists
            if fg_id in (i[0], i[1])
        ])
        short_distance = min(ds)
        short_distances.append(short_distance)

    return sum(short_distances)/len(short_distances)


def output_face_atoms(mol, metal_atom_ids, face_sets, prefix):
    """
    Output structure with metal atoms, differentiated by face.

    """

    _face_atom_types = {
        '001': 'H',
        '100': 'O',
        '010': 'F',
        '00-1': 'N',
        '-100': 'C',
        '0-10': 'B',
    }

    atom_count = 0
    xyz_file = f'{prefix}_metal_faces.xyz'
    atoms = []
    for fs in face_sets.sets:
        fsv = face_sets.vertices[fs]
        fs_atom_ids = tuple(metal_atom_ids[i] for i in fsv)
        fs_atom_positions = mol.get_atomic_positions(fs_atom_ids)

        for ai, ap in zip(fs_atom_ids, fs_atom_positions):
            atoms.append((
                _face_atom_types[fs],
                ap[0],
                ap[1],
                ap[2],
                ai,
            ))
            atom_count += 1

    with open(xyz_file, 'w') as f:
        f.write(f'{atom_count}\n\n')
        for ma in atoms:
            f.write(f'{ma[0]} {ma[1]} {ma[2]} {ma[3]} {ma[4]}\n')


def calculate_interior_face_angles(mol, metal_atom_ids, face_sets):
    """
    Calculate interor angles in faces of metal atoms.

    """

    pos_mat = mol.get_position_matrix()

    metal_metal_vectors = {
        (idx1, idx2): pos_mat[idx2]-pos_mat[idx1]
        for idx1, idx2 in permutations(metal_atom_ids, r=2)
    }

    cube_angles = {fs: {} for fs in face_sets.sets}
    for fs in face_sets.sets:
        fsv = face_sets.vertices[fs]
        fsc = face_sets.connected[fs]
        fs_atom_ids = tuple(metal_atom_ids[i] for i in fsv)

        # Calculate the interior angles based on connected metals.
        interior_angles = {}
        for idx in fsv:
            conn = [
                j
                for i in fsc if idx in i
                for j in i if j != idx
            ]
            fs_idx = metal_atom_ids[idx]
            fs_conn = [metal_atom_ids[i] for i in conn]

            # Define angle based on the two connections.
            pair1 = (fs_idx, fs_conn[0])
            pair2 = (fs_idx, fs_conn[1])
            # print(idx, pair1, pair2)
            vector1 = metal_metal_vectors[(pair1)]
            vector2 = metal_metal_vectors[(pair2)]
            interior_angle = np.degrees(
                angle_between(vector1, vector2)
            )
            if interior_angle < 75 or interior_angle > 105:
                print(
                    f'Interior angle: {interior_angle} - suggests bad '
                    'optimisation.'
                )
            interior_angles[idx] = interior_angle

        cube_angles[fs] = interior_angles

    return cube_angles


def convert_symm_names(symm_name=None, no_symbol=False):

    if no_symbol:
        new_names = {
            'd2': r'$D_\mathrm{2}$1',
            'th1': r'$T_\mathrm{h}$1',
            'th2': r'$T_\mathrm{h}$2',
            'td': r'$T$1',
            # 'tl': r'$T_\mathrm{\Lambda}$',
            's41': r'$S_\mathrm{4}$1',
            's42': r'$S_\mathrm{4}$2',
            's61': r'$S_\mathrm{6}$1',
            's62': r'$S_\mathrm{6}$2',
            'd31': r'$D_\mathrm{3}$1',
            'd32': r'$D_\mathrm{3}$2',
            'd31n': r'$D_\mathrm{3}$1n',
            'd32n': r'$D_\mathrm{3}$2n',
            'c2v': r'$C_\mathrm{2h}$1',
            'c2h': r'$C_\mathrm{2v}$1',
        }
    else:
        new_names = {
            'd2': r'$D_\mathrm{2}$',
            'th1': r'$T_\mathrm{h}$1',
            'th2': r'$T_\mathrm{h}$2',
            'td': r'$T$1-$\mathrm{\Delta}$',
            'tl': r'$T$1-$\mathrm{\Lambda}$',
            's41': r'$S_\mathrm{4}$1',
            's42': r'$S_\mathrm{4}$2',
            's61': r'$S_\mathrm{6}$1',
            's62': r'$S_\mathrm{6}$2',
            'd31': r'$D_\mathrm{3}$1',
            'd32': r'$D_\mathrm{3}$2',
            'd31n': r'$D_\mathrm{3}$1n',
            'd32n': r'$D_\mathrm{3}$2n',
            'c2v': r'$C_\mathrm{2h}$',
            'c2h': r'$C_\mathrm{2v}$',
        }
    if symm_name is None:
        return new_names
    else:
        return new_names[symm_name]


def convert_lig_names_from_cage(lig_name, as_int=False, as_sub=False):

    if as_int:
        new_names = {
            'quad2_5': 1,
            'quad2_16': 2,
            'quad2_12': 3,
            'quad2_3': 4,
            'quad2_8': 5,
            'quad2_2': 6,
            # 'quad2_1': 7,
        }
    elif as_sub:
        new_names = {
            'quad2_5': '$\\bf{A}$',  # '1',
            'quad2_16': '$\\bf{B}$',  # '2',
            'quad2_12': '$\\bf{C}$',  # '3',
            'quad2_3': '$\\bf{D}$',  # '4',
            'quad2_8': '$\\bf{E}$',  # '5',
            'quad2_2': '$\\bf{F}$',  # '6',
            # 'quad2_1': 7,
        }
    else:
        new_names = {
            'quad2_5': '$\\bf{A}$-Br',  # '1',
            'quad2_16': '$\\bf{B}$-Br',  # '2',
            'quad2_12': '$\\bf{C}$-Br',  # '3',
            'quad2_3': '$\\bf{D}$-Br',  # '4',
            'quad2_8': '$\\bf{E}$-Br',  # '5',
            'quad2_2': '$\\bf{F}$-Br',  # '6',
            # 'quad2_1': 'G',  # '7',
        }

    return new_names[lig_name]


def get_planar_conformer(molecule, N=100):
    cids, confs = build_conformers(
        mol=molecule,
        N=N,
        ETKDG_version='v3'
    )
    print('getting optimal conformer...')
    min_plane_dev = 100000000
    min_cid = -10

    new_molecule = molecule.clone()

    for cid in cids:

        # Update stk_mol to conformer geometry.
        new_molecule = update_from_rdkit_conf(
            stk_mol=new_molecule,
            rdk_mol=confs,
            conf_id=cid
        )

        plane_dev = calculate_molecule_planarity(new_molecule)
        if plane_dev < min_plane_dev:
            min_cid = cid
            min_plane_dev = plane_dev
            molecule = update_from_rdkit_conf(
                stk_mol=molecule,
                rdk_mol=confs,
                conf_id=min_cid
            )

    return molecule


def planarfy(ligands):
    """
    Get the most planar conformer of each ligand.

    This is done by determining the ETKDG conformer with the smallest
    plane deviation from its plane of best fit.

    """

    new_ligands = {}

    for ligand in ligands:
        planar_file = f'{ligand}_planar.mol'
        if os.path.exists(planar_file):
            opt_lig = ligands[ligand].with_structure_from_file(
                planar_file
            )
        else:
            print(f'doing {ligand}...')
            opt_lig = get_planar_conformer(ligands[ligand])
            opt_lig.write(planar_file)
        new_ligands[ligand] = opt_lig

    return new_ligands


def split_xyz_file(num_atoms, xyz_file):
    """
    Splits xyz trajectory file into xyz files.

    """

    with open(xyz_file, 'r') as f:
        lines = f.readlines()

    file_strings = []
    string = []
    for line in lines:
        if f' {num_atoms} ' in f' {line.strip()} ':
            if len(string) == 0:
                string.append(line)
            else:
                # New block.
                file_strings.append(string)
                string = [line]
        else:
            string.append(line)
    # Add last set.
    file_strings.append(string)

    out_files = []
    for i, fs in enumerate(file_strings):
        file_name = xyz_file.replace('.xyz', f'_s{i}.xyz')
        with open(file_name, 'w') as f:
            for line in fs:
                f.write(line)
        out_files.append(file_name)

    return out_files


def start_at_0(data_dict):

    new_dict = {}
    min_val = min([i for i in data_dict.values() if i is not None])
    for i in data_dict:
        if data_dict[i] is None:
            new_dict[i] = None
        else:
            new_dict[i] = data_dict[i] - min_val
    return new_dict


def get_plottables(measures, name):

    plottables = {
        'octop': {
            'data': measures['octop'],
            'ylabel': r'min. $q_{\mathrm{oct}}$',
            'ylim': (0, 1),
            'filename': f'{name}_minOPs.pdf'
        },
        'm_cube_shape': {
            'data': measures['m_cube_shape'],
            'ylabel': 'CU-8 shape measure',
            'ylim': (0, 1),
            'filename': f'{name}_cubeshape.pdf'
        },
        'lsesum': {
            'data': start_at_0(data_dict=measures['lsesum']),
            'ylabel': r'rel. sum strain energy [kJmol$^{-1}$]',
            'ylim': (0, 500),
            'filename': f'{name}_sumLSE.pdf'
        },
        'minitors': {
            'data': measures['minitors'],
            'ylabel': r'min. imine torsion [degrees]',
            'ylim': (0, 185),
            'filename': f'{name}_mintors.pdf'
        },
        'maxcrplan': {
            'data': measures['maxcrplan'],
            'ylabel': r'max. core planarity [$\mathrm{\AA}$]',
            'ylim': (0, 185),
            'filename': f'{name}_maxcrplane.pdf'
        },
        'maxintangledev': {
            'data': measures['maxintangledev'],
            'ylabel': r'max. interior angle deviation [$^{\circ}$]',
            'ylim': (-0.5, 10),
            'filename': f'{name}_maxintangledev.pdf'
        },
        'maxMLlength': {
            'data': measures['maxMLlength'],
            'ylabel': r'max. N-Zn bond length [$\mathrm{\AA}$]',
            'ylim': (2, 2.5),
            'filename': f'{name}_maxmld.pdf'
        },
        'porediam': {
            'data': measures['porediam'],
            'ylabel': r'pore diamater [$\mathrm{\AA}$]',
            'ylim': (0, 20),
            'filename': f'{name}_porediam.pdf'
        },
        'formatione': {
            'data': start_at_0(data_dict=measures['formatione']),
            'ylabel': r'rel. formation energy [kJmol$^{-1}$]',
            'ylim': (-10, 1000),
            'filename': f'{name}_relfe.pdf'
        },
    }

    return plottables


def calculate_cube_shape_measure(name, molecule):
    """
    Calculate the shape of an 8 atom molecule.

    Shape: http://www.ee.ub.edu/index.php?option=com_content&view=
    article&id=575:shape-available&catid=80:news&Itemid=466

    """

    if molecule.get_num_atoms() != 8:
        raise ValueError('Molecule does not have 8 atoms.')

    shape_dicts = (
        ref_shape_dict()['cube'],
        ref_shape_dict()['octagon']
    )
    n_verts = list(set([i['vertices'] for i in shape_dicts]))
    if len(n_verts) != 1:
        raise ValueError('Different vertex shapes selected.')

    input_file = f'{name}_shp.dat'
    std_out = f'{name}_shp.out'
    output_file = f'{name}_shp.tab'
    write_shape_input_file(
        input_file=input_file,
        name=name,
        structure=molecule,
        num_vertices=n_verts[0],
        central_atom_id=0,
        ref_shapes=[i['code'] for i in shape_dicts],
    )

    run_shape(input_file, env_set.shape_path(), std_out)
    shapes = collect_all_shape_values(output_file)
    return shapes


def ref_shape_dict():
    return {
        'cube': {
            'vertices': '8',
            'label': 'CU-8',
            'code': '4',
        },
        'octagon': {
            'vertices': '8',
            'label': 'OP-8',
            'code': '1',
        },
    }


def write_shape_input_file(
    input_file,
    name,
    structure,
    num_vertices,
    central_atom_id,
    ref_shapes,
):
    """
    Write input file for shape.

    """

    title = '$shape run by Andrew Tarzia.\n'
    size_of_poly = f'{num_vertices} {central_atom_id}\n'
    codes = ' '.join(ref_shapes)+'\n'

    structure_string = f'{name}\n'
    pos_mat = structure.get_position_matrix()
    for atom in structure.get_atoms():
        ele = atom.__class__.__name__
        x, y, z = pos_mat[atom.get_id()]
        structure_string += f'{ele} {x} {y} {z}\n'

    string = title+size_of_poly+codes+structure_string

    with open(input_file, 'w') as f:
        f.write(string)


def run_shape(input_file, shape_path, std_out):
    """
    Run input file for shape.

    """

    cmd = f'{shape_path} {input_file}'
    with open(std_out, 'w') as f:
        # Note that sp.call will hold the program until completion
        # of the calculation.
        sp.call(
            cmd,
            stdin=sp.PIPE,
            stdout=f,
            stderr=sp.PIPE,
            # Shell is required to run complex arguments.
            shell=True
        )


def collect_all_shape_values(output_file):
    """
    Collect shape values from output.

    """

    with open(output_file, 'r') as f:
        lines = f.readlines()

    label_idx_map = {}
    for line in reversed(lines):
        if 'Structure' in line:
            line = [
                i.strip()
                for i in line.rstrip().split(']')[1].split(' ')
                if i.strip()
            ]
            for idx, symb in enumerate(line):
                label_idx_map[symb] = idx
            break
        line = [i.strip() for i in line.rstrip().split(',')]
        values = line

    shapes = {
        i: float(values[1+label_idx_map[i]]) for i in label_idx_map
    }

    return shapes


def update_from_rdkit_conf(stk_mol, rdk_mol, conf_id):
    """
    Update the structure to match `conf_id` of `mol`.

    Parameters
    ----------
    struct : :class:`stk.Molecule`
        The molecule whoce coordinates are to be updated.

    mol : :class:`rdkit.Mol`
        The :mod:`rdkit` molecule to use for the structure update.

    conf_id : :class:`int`
        The conformer ID of the `mol` to update from.

    Returns
    -------
    :class:`.Molecule`
        The molecule.

    """

    pos_mat = rdk_mol.GetConformer(id=conf_id).GetPositions()
    return stk_mol.with_position_matrix(pos_mat)


def build_conformers(mol, N, ETKDG_version=None):
    """
    Convert stk mol into RDKit mol with N conformers.

    ETKDG_version allows the user to pick their choice of ETKDG params.

    `None` provides the settings used in ligand_combiner and unsymm.

    Other options:
        `v3`:
            New version from DOI: 10.1021/acs.jcim.0c00025
            with improved handling of macrocycles.

    """
    molecule = mol.to_rdkit_mol()
    molecule.RemoveAllConformers()

    if ETKDG_version is None:
        cids = rdkit.EmbedMultipleConfs(
            mol=molecule,
            numConfs=N,
            randomSeed=1000,
            useExpTorsionAnglePrefs=True,
            useBasicKnowledge=True,
            numThreads=4,
        )

    elif ETKDG_version == 'v3':
        params = rdkit.ETKDGv3()
        params.randomSeed = 1000
        cids = rdkit.EmbedMultipleConfs(
            mol=molecule,
            numConfs=N,
            params=params
        )

    print(f'there are {molecule.GetNumConformers()} conformers')
    return cids, molecule


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
        atom_ids=next(
            molecule.get_functional_groups()
        ).get_placer_ids(),
    )
    edge_position = target_coords[0]
    molecule = (
        molecule.with_rotation_to_minimize_angle(
            start=fg_bonder_centroid - centroid_pos,
            target=edge_position - edge_centroid,
            axis=edge_normal,
            origin=centroid_pos,
        )
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
    molecule = (
        molecule.with_rotation_to_minimize_angle(
            start=long_axis_vector,
            target=[1, 0, 0],
            axis=edge_normal,
            origin=centroid_pos,
        )
    )
    return molecule
