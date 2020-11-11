#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Modules defining and building the Cage class and subclasses.

Author: Andrew Tarzia

Date Created: 03 Mar 2020

"""

from copy import deepcopy
from os.path import exists
import json
import pywindow as pw
from os import system, mkdir

import stk

import atools
from molecule_building import (
    metal_FFs,
    optimize_SCA_complex,
    get_lowest_energy_conformer,
)
from face_sets import M8L6_FaceSets
from utilities import (
    calculate_paired_face_anisotropies,
    calculate_cube_likeness,
)


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
        self.bb_file = f'{self.name}_BBs'
        self.crush_file = f'{self.name}_cru'
        self.uff4mof_file = f'{self.name}_uff'
        self.uff4mof_CG_file = f'{self.name}_uffCG'
        self.uffMD_file = f'{self.name}_prextb'
        self.opt_file = f'{self.name}_optc'
        self.pw_file = f'{self.name}_pw'
        self.op_file = f'{self.name}_OP'
        self.cl_file = f'{self.name}_cl'
        self.ls_file = f'{self.name}_LSE'
        self.charge = charge
        self.free_electron_options = free_electron_options
        self.cage_set_dict = cage_set_dict

    def build(self):
        print(f'....building {self.name}')
        cage = stk.ConstructedMolecule(self.topology_graph)
        if not exists(f'{self.unopt_file}.mol'):
            cage.write(f'{self.unopt_file}.mol')
        self.cage = cage

    def save_bb_xyz(self):

        raise NotImplementedError('Currently broken.')

        if exists(f'{self.opt_file}.mol'):
            self.cage = self.cage.with_structure_from_file(
                f'{self.opt_file}.mol'
            )
        self.cage.write(f'{self.bb_file}.xyz')

        bb_pos = {
            stk.Smiles().get_key(i): j
            for i, j in self.building_blocks.items()
        }
        bb_data = {}
        for bb in self.cage.get_building_blocks():
            smiles = stk.Smiles().get_key(bb)
            if smiles not in bb_data:
                bb_data[smiles] = {
                    'no': len(bb_data),
                    'pos': bb_pos[smiles]
                }

        # Add column to XYZ file.
        with open(f'{self.bb_file}.xyz', 'r') as f:
            lines = f.readlines()

        new_lines = [i.rstrip() for i in lines]
        unique_ids = {}
        for i, nl in enumerate(new_lines):
            if i < 2:
                continue
            atom_id = i-2
            atom_info, = self.cage.get_atom_infos(atom_id)
            building_block = atom_info.get_building_block()
            building_block_id = atom_info.get_building_block_id()
            smi = stk.Smiles().get_key(building_block)
            bb_d = bb_data[smi]
            va = self.vertex_alignments[
                building_block_id
            ]
            bb_type = f"{bb_d['no']+1}{va+1}"
            if bb_type not in unique_ids:
                unique_ids[bb_type] = str(len(unique_ids)+1)

            bb_ids = unique_ids[bb_type]
            new_line = nl+' '+bb_ids
            new_lines[i] = new_line

        with open(f'{self.bb_file}.xyz', 'w') as f:
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
            return
        print(f'....optimizing {self.name}')

        # Run if crush output does not exist.
        if not exists(f'{self.crush_file}.mol'):
            self.cage = atools.MOC_collapse(
                cage=self.cage,
                cage_name=self.name,
                step_size=step_size,
                distance_cut=distance_cut,
                scale_steps=scale_steps,
            )
            self.cage.write(f'{self.crush_file}.mol')
        else:
            self.cage = self.cage.with_structure_from_file(
                f'{self.crush_file}.mol'
            )

        # Run if uff4mof opt output does not exist.
        if not exists(f'{self.uff4mof_CG_file}.mol'):
            self.cage = atools.MOC_uff_opt(
                self.cage,
                self.name,
                metal_FFs=custom_metal_FFs,
                CG=True,
                maxcyc=1000,
                metal_ligand_bond_order='',
            )
            self.cage.write(f'{self.uff4mof_CG_file}.mol')
        else:
            self.cage = self.cage.with_structure_from_file(
                f'{self.uff4mof_CG_file}.mol'
            )

        # Run if uff4mof opt output does not exist.
        if not exists(f'{self.uff4mof_file}.mol'):
            self.cage = atools.MOC_uff_opt(
                self.cage,
                self.name,
                metal_FFs=custom_metal_FFs,
                metal_ligand_bond_order='',
            )
            self.cage.write(f'{self.uff4mof_file}.mol')
        else:
            self.cage = self.cage.with_structure_from_file(
                f'{self.uff4mof_file}.mol'
            )

        # Run if uff4mof MD output does not exist.
        if not exists(f'{self.uffMD_file}.mol'):
            self.cage = atools.MOC_MD_opt(
                self.cage,
                self.name,
                integrator='leapfrog verlet',
                temperature=400,
                N=10,
                timestep=0.5,
                equib=0.1,
                production=2,
                metal_FFs=custom_metal_FFs,
                metal_ligand_bond_order='half',
                opt_conf=False,
                save_conf=False
            )
            self.cage.write(f'{self.uffMD_file}.mol')
        else:
            self.cage = self.cage.with_structure_from_file(
                f'{self.uffMD_file}.mol'
            )

        self.cage = atools.MOC_xtb_opt(
            self.cage,
            self.name,
            gfn_exec='/home/atarzia/software/xtb-6.3.1/bin/xtb',
            nc=6,
            free_e=free_e,
            charge=self.charge,
            opt_level='normal',
            etemp=300,
            solvent=(self.cage_set_dict['solvent'], 'normal'),
        )
        self.cage.write(f'{self.opt_file}.mol')

    def compare_UHF_values(self):
        print(f'....comparing UHF of {self.name}')
        raise NotImplementedError()

    def calculate_formation_energy(
        self,
        org_ligs,
        smiles_keys,
        file_prefix,
        gfn_exec,
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

        if self.topology_string != 'm8l6face':
            raise NotImplementedError('need to handle het cages.')

        components = deepcopy(self.cage_set_dict['components'])
        solvent = (self.cage_set_dict['solvent'], 'normal')

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
                        settings = {
                            'conf_opt_level': 'crude',
                            'final_opt_level': 'extreme',
                            'charge': charge,
                            'no_unpaired_e': no_unpaired_e,
                            'max_runs': 1,
                            'calc_hessian': False,
                            'solvent': solvent,
                            'crest_exec': (
                                '/home/atarzia/software/crest/crest'
                            ),
                            'nc': 4,
                            'etemp': 300,
                            'keepdir': False,
                            'cross': True,
                            'md_len': None,
                            'ewin': 5,
                            'speed_setting': 'squick',
                        }
                        conf_folder = (
                            f"{components[comp]['name']}_confs/"
                        )
                        if not exists(conf_folder):
                            mkdir(conf_folder)
                        temp_mol = get_lowest_energy_conformer(
                            name=components[comp]['name'],
                            mol=temp_mol,
                            settings=settings,
                            gfn_exec=(
                                '/home/atarzia/software/xtb-6.3.1/bin/'
                                'xtb'
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
                atools.calculate_energy(
                    name=f'{self.name}_{comp}',
                    mol=components[comp]['mol'],
                    gfn_exec=(
                        '/home/atarzia/software/xtb-6.3.1/bin/xtb'
                    ),
                    ey_file=ey_file,
                    charge=charge,
                    no_unpaired_e=no_unpaired_e,
                    solvent=solvent
                )
            ey = atools.read_gfnx2xtb_eyfile(ey_file)
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
    ):
        """
        Analyse cage ligand geometry for strain.

        """

        print(f'....analyzing ligand geometry of {self.name}')
        # Collect the atomic positions of the organic linkers in the
        # cage for analysis.
        org_ligs, smiles_keys = atools.get_organic_linkers(
            cage=self.cage,
            metal_atom_nos=(metal_atom_no, ),
            file_prefix=f'{self.name}_sg'
        )

        num_unique_ligands = len(set(smiles_keys.values()))
        if num_unique_ligands != expected_ligands:
            raise UnexpectedNumLigands(
                f'{self.name} had {num_unique_ligands} unique ligands'
                f', {expected_ligands} were expected. Suggests bad '
                'optimization.'
            )

        atools.get_lowest_energy_conformers(
            org_ligs=org_ligs,
            smiles_keys=smiles_keys,
            file_prefix=f'{self.base_name}_sg',
            gfn_exec='/home/atarzia/software/xtb-6.3.1/bin/xtb',
            conformer_function=get_lowest_energy_conformer,
            conformer_settings={
                'conf_opt_level': 'crude',
                'final_opt_level': 'extreme',
                'charge': 0,
                'no_unpaired_e': 0,
                'max_runs': 1,
                'calc_hessian': False,
                'solvent': (self.cage_set_dict['solvent'], 'normal'),
                'crest_exec': '/home/atarzia/software/crest/crest',
                'nc': 4,
                'etemp': 300,
                'keepdir': False,
                'cross': True,
                'md_len': None,
                'ewin': 5,
                'speed_setting': 'squick',
            },
        )

        self.calculate_formation_energy(
            org_ligs=org_ligs,
            smiles_keys=smiles_keys,
            file_prefix=f'{self.base_name}_sg',
            gfn_exec='/home/atarzia/software/xtb-6.3.1/bin/xtb',
            cage_free_e=free_e,
        )

        lse_dict = atools.calculate_ligand_SE(
            org_ligs=org_ligs,
            smiles_keys=smiles_keys,
            output_json=f'{self.ls_file}.json',
            file_prefix=f'{self.base_name}_sg'
        )
        imine_torsion_dict = atools.calculate_abs_imine_torsions(
            org_ligs=org_ligs
        )
        planarity_dict = atools.calculate_ligand_planarities(
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
                op_set = atools.get_order_values(
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
        self.bl_data = atools.calculate_metal_ligand_distance(
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
                pass

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
