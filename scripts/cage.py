#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Modules defining and building the Cage class and subclasses.

Author: Andrew Tarzia

Date Created: 03 Mar 2020

"""

from copy import deepcopy
from os.path import exists, join
import numpy as np
import json
from itertools import combinations
import pywindow as pw
from os import system, mkdir

import stk
import stko

from molecule_building import (
    metal_FFs,
    optimize_SCA_complex,
    get_lowest_energy_conformer,
)
from face_sets import M8L6_FaceSets
from utilities import (
    calculate_paired_face_anisotropies,
    calculate_cube_likeness,
    calculate_cube_shape_measure,
    calculate_abs_imine_torsions,
    calculate_ligand_SE,
    calculate_metal_ligand_distance,
    calculate_ligand_planarities,
    get_query_atom_ids,
    get_order_values,
    calculate_energy,
    read_gfnx2xtb_eyfile,
    get_organic_linkers,
    get_lowest_energy_conformers,
)
import env_set


class UnexpectedNumLigands(Exception):
    ...


class PrecursorNotOptimizedError(Exception):
    ...


class CageNotOptimizedError(Exception):
    ...


class Cage:
    """
    Generic class that builds and analyses stk.ConstructuedMolecules.

    """

    def __init__(
        self,
        name,
        base_name,
        topology_fn,
        topology_string,
        symmetry_string,
        building_blocks,
        vertex_alignments,
        charge,
        free_electron_options,
        cage_set_dict
    ):

        self.name = name
        self.base_name = base_name
        self.topology_fn = topology_fn
        self.topology_string = topology_string
        self.symmetry_string = symmetry_string
        self.building_blocks = building_blocks
        self.vertex_alignments = vertex_alignments
        self.topology_graph = self.topology_fn(
            building_blocks=self.building_blocks,
            vertex_alignments=self.vertex_alignments,
        )
        self.unopt_file = f'{self.name}_unopt'
        self.crush_file = f'{self.name}_cru'
        self.uff4mof_file = f'{self.name}_uff'
        self.uff4mof_CG_file = f'{self.name}_uffCG'
        self.uffMD_file = f'{self.name}_prextb'
        self.opt_file = f'{self.name}_optc'
        self.pw_file = f'{self.name}_pw'
        self.op_file = f'{self.name}_OP'
        self.cl_file = f'{self.name}_cl'
        self.sh_file = f'{self.name}_sh'
        self.ls_file = f'{self.name}_LSE'
        self.charge = charge
        self.free_electron_options = free_electron_options
        self.cage_set_dict = cage_set_dict
        self.optimized = None

    def build(self):
        print(f'....building {self.name}')
        cage = stk.ConstructedMolecule(self.topology_graph)
        if not exists(f'{self.unopt_file}.mol'):
            cage.write(f'{self.unopt_file}.mol')
        self.cage = cage

    def save_bb_vector_xyzs(self, structure_file):
        """
        Save an XYZ for visualisation of the cage at an opt stage.

        """

        def get_ligand_vectors(mol, metal_atom_nos):
            """
            Define ligand vectors based on bonder positions.

            """

            l_vectors = []

            # Get building block ids that do not have metal atoms.
            bb_atom_ids = {}
            metal_bbs = []
            for ainfo in mol.get_atom_infos():
                bbid = ainfo.get_building_block_id()
                if bbid not in bb_atom_ids:
                    bb_atom_ids[bbid] = []
                bb_atom_ids[bbid].append(ainfo.get_atom().get_id())
                a_atom_nu = ainfo.get_atom().get_atomic_number()
                if a_atom_nu in metal_atom_nos:
                    metal_bbs.append(bbid)

            ligand_bb_atom_ids = {
                i: bb_atom_ids[i]
                for i in bb_atom_ids if i not in metal_bbs
            }

            # Add centroid to l_vectors.
            for bb in ligand_bb_atom_ids:
                x, y, z = mol.get_centroid(ligand_bb_atom_ids[bb])
                l_vectors.append([[x, y, z, 'He']])

            # Get axis vectors.
            # Use Smarts to find connections between building blocks.
            # Which is used to find bonder atoms on ligands.
            # NX3 - CX3, where X != H.
            all_lig_ids = set([
                j for i in ligand_bb_atom_ids.values() for j in i
            ])
            smarts = '[#7X3]~[#6X3]'
            rdkit_mol = mol.to_rdkit_mol()
            query_ids = get_query_atom_ids(smarts, rdkit_mol)
            bonder_atom_ids = []
            for atom_ids in query_ids:
                n_id = atom_ids[0]
                c_id = atom_ids[1]
                if n_id not in all_lig_ids and c_id in all_lig_ids:
                    bonder_atom_ids.append(c_id)

            for bb in ligand_bb_atom_ids:
                bb_bonders = [
                    i for i in bonder_atom_ids
                    if i in ligand_bb_atom_ids[bb]
                ]
                for atom_ids in combinations(bb_bonders, r=2):
                    c1_pos = tuple(
                        mol.get_atomic_positions(atom_ids[0])
                    )[0]
                    c2_pos = tuple(
                        mol.get_atomic_positions(atom_ids[1])
                    )[0]
                    vector = c2_pos - c1_pos
                    pts = [
                        np.linspace(c2_pos[i], c1_pos[i], 10)
                        for i in np.arange(len(c2_pos))
                    ]
                    vector = [[i, j, k, 'Kr'] for i, j, k in zip(*pts)]
                    l_vectors.append(vector)

            return l_vectors

        def get_metal_vectors(mol, metal_atom_nos):
            """
            Define metal vectors based on N-N vectors in each complex.

            """

            m_vectors = []
            pos_mat = mol.get_position_matrix()

            # Add a single point vector for each metal atom.
            for atom in mol.get_atoms():
                if atom.get_atomic_number() in metal_atom_nos:
                    x, y, z = pos_mat[atom.get_id()]
                    m_vectors.append([[x, y, z, 'Ar']])

            # Add vectors for all shortest N-N lengths.
            # Use Smarts to find pairs of N atoms.
            # N(X3)-CX2H1-CX3-NX3, where X != H.
            smarts = '[#7X3]~[#6X3H1]~[#6X3!H1]~[#7X3]'
            rdkit_mol = mol.to_rdkit_mol()
            query_ids = get_query_atom_ids(smarts, rdkit_mol)
            for atom_ids in query_ids:
                n1_pos = tuple(
                    mol.get_atomic_positions(atom_ids[0])
                )[0]
                n2_pos = tuple(
                    mol.get_atomic_positions(atom_ids[3])
                )[0]
                pts = [
                    np.linspace(n2_pos[i], n1_pos[i], 5)
                    for i in np.arange(len(n2_pos))
                ]
                vector = [[i, j, k, 'Xe'] for i, j, k in zip(*pts)]
                m_vectors.append(vector)

            return m_vectors

        ve_output_file = structure_file.replace('.mol', '_VECTs.xyz')

        temp_mol = self.cage.with_structure_from_file(structure_file)

        # Write vectors out to xyz file.
        temp_mol.write(ve_output_file)
        # Do so by adding rows to XYZ file.
        # Define the sets of vectors.
        metal_atom_nos = list(set([
            i.get_atomic_number() for i in self.cage.get_atoms()
            if i.get_atomic_number() in metal_FFs(CN=4).keys()
        ]))
        ligand_vectors = get_ligand_vectors(temp_mol, metal_atom_nos)
        metal_vectors = get_metal_vectors(temp_mol, metal_atom_nos)
        with open(ve_output_file, 'r') as f:
            lines = f.readlines()
        new_lines = [i.rstrip() for i in lines]
        for lv in ligand_vectors:
            for point in lv:
                new_lines.append(
                    f'{point[3]} {point[0]} {point[1]} {point[2]}'
                )
        for mv in metal_vectors:
            for point in mv:
                new_lines.append(
                    f'{point[3]} {point[0]} {point[1]} {point[2]}'
                )

        new_lines[0] = str(len(new_lines)-2)
        with open(ve_output_file, 'w') as f:
            for line in new_lines:
                f.write(line+'\n')

    def optimize(
        self,
        free_e,
        step_size,
        distance_cut,
        scale_steps
    ):
        custom_metal_FFs = metal_FFs(CN=6)

        # Skip if _opt.mol exists.
        if exists(f'{self.opt_file}.mol'):
            self.cage = self.cage.with_structure_from_file(
                f'{self.opt_file}.mol'
            )
            self.optimized = True
            return
        print(f'....optimizing {self.name}')
        self.optimized = None

        # Run if crush output does not exist.
        if not exists(f'{self.crush_file}.mol'):
            print(f'..doing collapser optimisation of {self.name}')
            output_dir = f'cage_opt_{self.name}_coll'
            optimizer = stko.Collapser(
                output_dir=output_dir,
                step_size=step_size,
                distance_cut=distance_cut,
                scale_steps=scale_steps,
            )
            self.cage = optimizer.optimize(mol=self.cage)
            self.cage.write(f'{self.crush_file}.mol')
        else:
            self.cage = self.cage.with_structure_from_file(
                f'{self.crush_file}.mol'
            )

        # Run if uff4mof opt output does not exist.
        if not exists(f'{self.uff4mof_CG_file}.mol'):
            CG = True
            maxcyc = 1000
            metal_ligand_bond_order = ''
            output_dir = (
                f'cage_opt_{self.name}_uff' if CG is False
                else f'cage_opt_{self.name}_uffCG'
            )
            print(f'..doing UFF4MOF optimisation of {self.name}')
            print(f'Conjugate Gradient: {CG}, Max steps: {maxcyc}')
            gulp_opt = stko.GulpUFFOptimizer(
                gulp_path=env_set.gulp_path(),
                maxcyc=maxcyc,
                metal_FF=custom_metal_FFs,
                metal_ligand_bond_order=metal_ligand_bond_order,
                output_dir=output_dir,
                conjugate_gradient=CG
            )
            gulp_opt.assign_FF(self.cage)
            self.cage = gulp_opt.optimize(mol=self.cage)
            self.cage.write(f'{self.uff4mof_CG_file}.mol')
        else:
            self.cage = self.cage.with_structure_from_file(
                f'{self.uff4mof_CG_file}.mol'
            )

        # Run if uff4mof opt output does not exist.
        if not exists(f'{self.uff4mof_file}.mol'):
            CG = False
            maxcyc = 1000
            metal_ligand_bond_order = ''
            output_dir = (
                f'cage_opt_{self.name}_uff' if CG is False
                else f'cage_opt_{self.name}_uffCG'
            )
            print(f'..doing UFF4MOF optimisation of {self.name}')
            print(f'Conjugate Gradient: {CG}, Max steps: {maxcyc}')
            gulp_opt = stko.GulpUFFOptimizer(
                gulp_path=env_set.gulp_path(),
                maxcyc=maxcyc,
                metal_FF=custom_metal_FFs,
                metal_ligand_bond_order=metal_ligand_bond_order,
                output_dir=output_dir,
                conjugate_gradient=CG
            )
            gulp_opt.assign_FF(self.cage)
            self.cage = gulp_opt.optimize(mol=self.cage)
            self.cage.write(f'{self.uff4mof_file}.mol')
        else:
            self.cage = self.cage.with_structure_from_file(
                f'{self.uff4mof_file}.mol'
            )

        # Run if uff4mof MD output does not exist.
        if not exists(f'{self.uffMD_file}.mol'):
            print(f'..doing UFF4MOF MD of {self.name}')
            gulp_MD = stko.GulpUFFMDOptimizer(
                gulp_path=env_set.gulp_path(),
                metal_FF=custom_metal_FFs,
                metal_ligand_bond_order='half',
                output_dir=f'cage_opt_{self.name}_MD',
                integrator='leapfrog verlet',
                ensemble='nvt',
                temperature=400,
                equilbration=0.1,
                production=2,
                timestep=0.5,
                N_conformers=10,
                opt_conformers=False,
                save_conformers=False,
            )
            gulp_MD.assign_FF(self.cage)
            self.cage = gulp_MD.optimize(self.cage)
            self.cage.write(f'{self.uffMD_file}.mol')
        else:
            self.cage = self.cage.with_structure_from_file(
                f'{self.uffMD_file}.mol'
            )

        try:
            print(f'..........doing XTB optimisation of {self.name}')
            xtb_opt = stko.XTB(
                xtb_path=env_set.xtb_path(),
                output_dir=f'cage_opt_{self.name}_xtb',
                gfn_version=2,
                num_cores=6,
                opt_level='normal',
                charge=self.charge,
                num_unpaired_electrons=free_e,
                max_runs=1,
                electronic_temperature=300,
                calculate_hessian=False,
                unlimited_memory=True,
                solvent=self.cage_set_dict['solvent'],
            )
            self.cage = xtb_opt.optimize(mol=self.cage)
            self.cage.write(f'{self.opt_file}.mol')
            self.optimized = True
        except (stko.XTBConvergenceError, stko.XTBOptimizerError):
            # Check if the optimisation was even attempted.
            opt_output_file = join(
                f'cage_opt_{self.name}_xtb',
                'optimization_1.output'
            )
            if not exists(opt_output_file):
                # If not, raise error and exit.
                raise CageNotOptimizedError(
                    'xTB optimisation of cage not even attempted for'
                    f'{self.name}. Try rerunning and check xTB is '
                    'installed.'
                )

            # Check if that output file actually contains some
            # steps because xtb may fail on initialisation.
            steps_lines = []
            has_steps = False
            with open(opt_output_file, 'r') as f:
                for line in f.readlines():
                    if ' CYCLE ' in line:
                        steps_lines.append(line)
            if len(steps_lines) > 0:
                has_steps = True
            if has_steps:
                # Set optimized to False, this avoids all analysis.
                self.optimized = False
            else:
                # If not, raise error and exit.
                raise CageNotOptimizedError(
                    'xTB optimisation of cage not even attempted '
                    f'for {self.name}. Try rerunning and check xTB'
                    ' is installed.'
                )

    def compare_UHF_values(self):
        print(f'....comparing UHF of {self.name}')
        raise NotImplementedError()

    def calculate_formation_energy(
        self,
        org_ligs,
        smiles_keys,
        file_prefix,
        xtb_path,
        cage_free_e,
    ):
        """
        Calculate cage formation energy.

        Defined as (for a homoleptic system):
        FE = [(cage energy) + a*(free precursor lig. energy)]
             - [b*(metal precursor energy) + c*(free ligand energy)]

        Where a = 6*b (octahedral metal complex).

        In practice, this is done alchemically based on the building
        blocks present in the molecules.

        FE =
        [(cage energy) + a*(energy of deleter atoms)] -
        [sum_i(b_i*(energy BB_i))]

        Where a is determined by the atoms in each BB and i is over all
        BBs used to build the cage.

        """

        if self.topology_string not in ['m8l6face', 'm8l6knot']:
            raise NotImplementedError('need to handle het cages.')

        components = deepcopy(self.cage_set_dict['components'])
        solvent = self.cage_set_dict['solvent']

        # Get lowest energy conformer filenames.
        low_e_lig_filenames = []
        for lig in org_ligs:
            stk_lig = org_ligs[lig]
            smiles_key = stk.Smiles().get_key(stk_lig)
            idx = smiles_keys[smiles_key]
            sgt = str(stk_lig.get_num_atoms())
            # Get optimized ligand name that excludes any cage
            # information.
            if file_prefix is None:
                filename_ = f'organic_linker_s{sgt}_{idx}_opt.mol'
            else:
                filename_ = f'{file_prefix}{sgt}_{idx}_opt.mol'
            low_e_lig_filenames.append(filename_)
        low_e_lig_filenames = set(low_e_lig_filenames)
        print(low_e_lig_filenames)

        # Load in each component.
        for comp in components:
            print(comp, components[comp])
            if comp in ['mprec', 'mpreclig']:
                opt_file = f"{components[comp]['name']}_opt.mol"
                low_e_file = f"{components[comp]['name']}_loweopt.mol"
                charge = components[comp]['charge']
                no_unpaired_e = components[comp]['unpaired_e']
                ey_file = f"{components[comp]['name']}_opt.ey"
                # Get lowest energy conformer of mprec or mpreclig.
                if components[comp]['smiles'] is not None:
                    temp_mol = stk.BuildingBlock(
                        smiles=components[comp]['smiles'],
                    )
                    # Optimisation and lowest energy conformer search.
                    if exists(opt_file):
                        temp_mol = temp_mol.with_structure_from_file(
                            opt_file
                        )
                    else:
                        temp_mol = optimize_SCA_complex(
                            complex=temp_mol,
                            name=components[comp]['name'],
                            dict={
                                'total_charge': charge,
                                'unpaired_e': no_unpaired_e
                            },
                            metal_FFs=metal_FFs(CN=6)
                        )
                        temp_mol.write(opt_file)
                    if exists(low_e_file):
                        temp_mol = temp_mol.with_structure_from_file(
                            low_e_file
                        )
                    else:
                        conf_folder = (
                            f"{components[comp]['name']}_confs/"
                        )
                        if not exists(conf_folder):
                            mkdir(conf_folder)
                        temp_mol = get_lowest_energy_conformer(
                            name=components[comp]['name'],
                            mol=temp_mol,
                            settings=env_set.crest_conformer_settings(
                                solvent=None
                            ),
                        )
                        temp_mol.write(low_e_file)

                elif exists(low_e_file):
                    temp_mol = stk.BuildingBlock.init_from_file(
                        low_e_file
                    )
                else:
                    raise PrecursorNotOptimizedError(
                        f'{low_e_file} not found!'
                    )
                components[comp]['mol'] = temp_mol

            elif comp in ['tritopic', 'tetratopic']:
                charge = components[comp]['charge']
                no_unpaired_e = components[comp]['unpaired_e']
                lig_name = self.cage_set_dict[comp]
                low_e_files = [
                    i for i in low_e_lig_filenames if lig_name in i
                ]
                if len(low_e_files) > 1:
                    raise UnexpectedNumLigands(
                        f'Found {len(low_e_files)} low energy '
                        f'conformers with name {lig_name}.'
                    )
                low_e_file = low_e_files[0]
                ey_file = low_e_file.replace('.mol', '.ey')
                temp_mol = stk.BuildingBlock.init_from_file(low_e_file)
                components[comp]['mol'] = temp_mol

            elif comp == 'cage':
                charge = self.charge
                no_unpaired_e = cage_free_e
                ey_file = f'{self.opt_file}.ey'
                if exists(f'{self.opt_file}.mol'):
                    components[comp]['mol'] = self.cage
                else:
                    raise CageNotOptimizedError(
                        'Expected cage to be optimized at this point: '
                        f'{self.opt_file}.mol does not exist.'
                    )

            # Calculate all components energies.
            if not exists(ey_file):
                calculate_energy(
                    name=f'{self.name}_{comp}',
                    mol=components[comp]['mol'],
                    xtb_path=env_set.xtb_path(),
                    ey_file=ey_file,
                    charge=charge,
                    no_unpaired_e=no_unpaired_e,
                    solvent=solvent
                )
            ey = read_gfnx2xtb_eyfile(ey_file)
            components[comp]['total_e'] = ey

        # Calculate formation energy in kJ/mol.
        prod_ey = sum([
            components[i]['total_e']*components[i]['count']
            for i in components if components[i]['product'] is True
        ])
        react_ey = sum([
            components[i]['total_e']*components[i]['count']
            for i in components if components[i]['product'] is False
        ])
        fe = prod_ey - react_ey
        print(prod_ey, react_ey, fe)
        self.fe_data = fe
        print(self.name, self.fe_data)

    def analyze_ligand_strain(
        self,
        metal_atom_no,
        expected_ligands,
        free_e,
        ligand_dir,
    ):
        """
        Analyse cage ligand geometry for strain.

        """

        print(f'....analyzing ligand geometry of {self.name}')
        # Collect the atomic positions of the organic linkers in the
        # cage for analysis.
        org_ligs, smiles_keys = get_organic_linkers(
            cage=self.cage,
            metal_atom_nos=(metal_atom_no, ),
            file_prefix=f'{self.name}_sg'
        )

        num_unique_ligands = len(set(smiles_keys.values()))
        if num_unique_ligands != expected_ligands:
            raise UnexpectedNumLigands(
                f'{self.name} had {num_unique_ligands} unique ligands'
                f', {expected_ligands} were expected. Suggests bad '
                'optimization. Recommend reoptimising structure.'
            )

        # Loads from flexibility analysis and optimizes with solvent.
        get_lowest_energy_conformers(
            org_ligs=org_ligs,
            smiles_keys=smiles_keys,
            file_prefix=f'{self.base_name}_sg',
            ligand_dir=ligand_dir,
            settings=env_set.crest_conformer_settings(
                solvent=self.cage_set_dict['solvent'],
            ),
        )

        self.calculate_formation_energy(
            org_ligs=org_ligs,
            smiles_keys=smiles_keys,
            file_prefix=f'{self.base_name}_sg',
            xtb_path=env_set.xtb_path,
            cage_free_e=free_e,
        )

        lse_dict = calculate_ligand_SE(
            org_ligs=org_ligs,
            smiles_keys=smiles_keys,
            output_json=f'{self.ls_file}.json',
            file_prefix=f'{self.base_name}_sg'
        )
        imine_torsion_dict = calculate_abs_imine_torsions(
            org_ligs=org_ligs,
        )
        for ol in imine_torsion_dict:
            if len(imine_torsion_dict[ol]) != 4:
                raise ValueError(
                  f'{len(imine_torsion_dict[ol])} minies found, '
                  'but 4 expected.'
                )
        planarity_dict = calculate_ligand_planarities(
            org_ligs=org_ligs
        )

        self.ls_data = {
            'strain_energies': lse_dict,
            'imine_torsions': imine_torsion_dict,
            'core_planarities': planarity_dict
        }

    def get_cage_face_sets(self):

        if self.topology_string == 'm8l6face':
            return M8L6_FaceSets(self.symmetry_string)
        else:
            return None

    def analyze_cube_shape(self):
        """
        Analyse cage geometry based on its `cube-likeness`.

        Cube likeness is defined by the metal positions.

        """

        if self.topology_string != 'm8l6face':
            raise NotImplementedError(
                f'Cube-likeness is not defined for '
                f'{self.topology_string} topology.'
            )

        print(f'....analyzing cube shape of {self.name}')

        # Calculate cube shape measure.
        # Metal atom structure.
        # Get metal atom ids.
        metal_atom_ids = [
            i.get_id() for i in self.cage.get_atoms()
            if i.get_atomic_number() in metal_FFs(CN=4).keys()
        ]
        # Write to mol file.
        self.cage.write(f'{self.name}_M.mol', atom_ids=metal_atom_ids)
        m_structure = stk.BuildingBlock.init_from_file(
            f'{self.name}_M.mol'
        )
        shapes = calculate_cube_shape_measure(self.name, m_structure)
        return shapes['CU-8']

    def analyze_cube_likeness(self):
        """
        Analyse cage geometry based on its `cube-likeness`.

        Cube likeness is defined by the metal positions.

        """

        if self.topology_string != 'm8l6face':
            raise NotImplementedError(
                f'Cube-likeness is not defined for '
                f'{self.topology_string} topology.'
            )

        # Check if output file exists.
        if not exists(f'{self.cl_file}.json'):
            print(f'....analyzing cube likeness of {self.name}')

            # Get metal-metal distances and face anisotropies.
            # Assumes that atom ordering of metals follows vertex
            # ordering in topology definition.

            # Get metal atom ids.
            metal_atom_ids = [
                i.get_id() for i in self.cage.get_atoms()
                if i.get_atomic_number() in metal_FFs(CN=4).keys()
            ]

            face_sets = self.get_cage_face_sets()
            self.fa_data = calculate_paired_face_anisotropies(
                mol=self.cage,
                metal_atom_ids=metal_atom_ids,
                face_sets=face_sets,
            )

            self.cl_data = calculate_cube_likeness(
                mol=self.cage,
                metal_atom_ids=metal_atom_ids,
                face_sets=face_sets,
            )

    def analyze_metal_strain(self):
        """
        Analyse cage geometry using order parameters.

        """

        # Check if output file exists.
        if not exists(f'{self.op_file}.json'):
            print(f'....analyzing metal geometry of {self.name}')

            # Get atomic numbers of all present metals.
            pres_atm_no = list(set([
                i.get_atomic_number() for i in self.cage.get_atoms()
                if i.get_atomic_number() in metal_FFs(CN=4).keys()
            ]))

            # Get OPs for each metal independantly.
            # Save to dict with atom id and atom type.
            op_res = {}
            for metal in pres_atm_no:
                op_set = get_order_values(
                    mol=self.cage,
                    metal=metal,
                    per_site=True
                )
                op_res[metal] = op_set

            # Save to output files.
            with open(f'{self.op_file}.json', 'w') as f:
                json.dump(op_res, f)

        # Get data.
        with open(f'{self.op_file}.json', 'r') as f:
            self.op_data = json.load(f)

        # Get metal-ligand binder atom bond length.
        print(
            'Warning! Calculate metal-ligand distance is hard-coded '
            'for Zn-N'
        )
        self.bl_data = calculate_metal_ligand_distance(
            mol=self.cage,
            metal_atomic_number=30,
            ligand_atomic_number=7,
        )

    def analyze_porosity(self, dump_molecule=False):
        """
        Analyse cage porosity with pywindow.

        """

        # Check if output file exists.
        if not exists(f'{self.pw_file}.json'):
            print(f'....analyzing porosity of {self.name}')
            # Load cage into pywindow.
            self.cage.write('temp.xyz')
            pw_cage = pw.MolecularSystem.load_file(
                'temp.xyz'
            )
            pw_cage_mol = pw_cage.system_to_molecule()
            system('rm temp.xyz')

            # Calculate pore size.
            try:
                pw_cage_mol.calculate_pore_diameter_opt()
                pw_cage_mol.calculate_pore_volume_opt()
            except ValueError:
                # Handle failure.
                pw_cage_mol.properties['pore_volume_opt'] = 0
                pw_cage_mol.properties['pore_diameter_opt'] = {
                    'diameter': 0,
                    'atom_1': 0,
                    'centre_of_mass': [0, 0, 0],
                }

            # Save files.
            pw_cage_mol.dump_properties_json(f'{self.pw_file}.json')
            if dump_molecule:
                pw_cage_mol.dump_molecule(
                    f'{self.pw_file}.pdb',
                    include_coms=True
                )

        # Get data.
        with open(f'{self.pw_file}.json', 'r') as f:
            self.pw_data = json.load(f)

    def __str__(self):
        return (
            f'{self.__class__.__name__}'
            f'(name={self.name}, topology={self.topology_graph})'
        )

    def __repr__(self):
        return str(self)
