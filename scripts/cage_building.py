#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Modules defining and building the Cage class and subclasses.

Author: Andrew Tarzia

Date Created: 03 Mar 2020

"""

from itertools import product
from os.path import exists, join
import matplotlib.pyplot as plt
import json
import pywindow as pw
import os

import stk

import atools
from molecule_building import metal_FFs
import symmetries
from utilities import calculate_binding_AR


def m4l4spacer_graph(metal_corner, ligand, ligand2=None):

    raise NotImplementedError('see m8l6 for example')

    if ligand.get_num_functional_groups() != 3:
        raise ValueError(f'{ligand} does not have 3 functional groups')

    cage = stk.ConstructedMolecule(
        topology_graph=stk.cage.M4L4Tetrahedron(
            building_blocks={
                metal_corner: range(4),
                ligand: range(4, 8)
            },
        )
    )

    return cage


def m8l6_graph(building_blocks, vertex_alignments):

    topology_graph = stk.cage.M8L6Cube(
        building_blocks=building_blocks,
        vertex_alignments=vertex_alignments,
    )

    return topology_graph


def m6l2l3_graph(metal_corner, ligand, ligand2):

    raise NotImplementedError('see m8l6 for example')

    lig_num_fgs = ligand.get_num_functional_groups()
    lig2_num_fgs = ligand2.get_num_functional_groups()

    if lig_num_fgs == 3 and lig2_num_fgs == 4:
        tri_ligand = ligand
        tet_ligand = ligand2
    elif lig_num_fgs == 4 and lig2_num_fgs == 3:
        tri_ligand = ligand2
        tet_ligand = ligand
    else:
        raise ValueError(
            f'{ligand} or {ligand2} does not have 3 or 4 functional '
            'groups'
        )

    cage = stk.ConstructedMolecule(
        topology_graph=stk.cage.M6L2L3Prism(
            building_blocks={
                metal_corner: range(6),
                tri_ligand: range(6, 8),
                tet_ligand: range(8, 11)
            },
        )
    )

    return cage


def available_topologies(string):
    """
    Get stk function of desired topology.

    """

    topologies = {
        'm4l4spacer': m4l4spacer_graph,
        'm8l6face': m8l6_graph,
        'm6l2l3': m6l2l3_graph,
    }

    try:
        return topologies[string]
    except KeyError:
        raise KeyError(f'{string} not in {topologies.keys()}')


class Cage:
    """
    Generic class that builds and analyses stk.ConstructuedMolecules.

    """

    def __init__(
        self,
        name,
        topology_fn,
        topology_string,
        building_blocks,
        vertex_alignments,
        charge,
        free_electron_options
    ):

        self.name = name
        self.topology_fn = topology_fn
        self.topology_string = topology_string
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
        self.ls_file = f'{self.name}_LSE'
        self.charge = charge
        self.free_electron_options = free_electron_options

    def build(self):
        print(f'....building {self.name}')
        cage = stk.ConstructedMolecule(self.topology_graph)
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

    def optimize(self, free_e, step_size, distance_cut, scale_steps):
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
                self.cage,
                self.name,
                step_size=step_size,
                distance_cut=distance_cut,
                scale_steps=scale_steps
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
                CG=True
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
                metal_FFs=custom_metal_FFs
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
                temperature=700,
                N=10,
                timestep=0.5,
                equib=0.1,
                production=2,
                metal_FFs=custom_metal_FFs,
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
            # solvent=('dmso', 'verytight')
        )
        self.cage.write(f'{self.opt_file}.mol')

    def compare_UHF_values(self):
        print(f'....comparing UHF of {self.name}')
        raise NotImplementedError()

    def analyze_formation_energy(self):
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
        raise NotImplementedError()

        # Define reactant list (all building blocks * count in cage)
        reactant_ey_files = []
        reactant_molecules = []
        # Iterate over all BBs and add to product molecules.
        # Define charge and no. unpaired e from lib file.

        print(reactant_ey_files, reactant_molecules)

        # Define product list (cage + deleters).
        cage = self.cage.update_from_file(f'{self.opt_file}.mol')
        product_ey_files = []
        product_molecules = [
            (cage, cage_name, cage_ey_file, cage_charge, cage_no_e)
        ]

        # Iterate over deleters and add to product molecules.
        # Define charge apriori, Br is -1.

        print(product_molecules, product_ey_files)

        # Run calculations with xTB.

        # Calculate FE in KJ/mol.
        react_eys = [read_ey(f'{i}.ey') for i in reactant_ey_files]
        produ_eys = [read_ey(f'{i}.ey') for i in product_ey_files]
        #  kJ/mol
        FE = produ_eys - react_eys
        self.FE = FE

    def get_organic_linkers(self, metal_atom_no):
        """
        Extract a list of organic linker structures from the cage.

        """

        org_lig = {}

        # Produce a graph from the cage that does not include metals.
        cage_g = nx.Graph()
        atom_ids_in_G = set()
        for atom in self.cage.atoms:
            if atom.atomic_number == metal_atom_no:
                continue
            cage_g.add_node(atom)
            atom_ids_in_G.add(atom.id)

        # Add edges.
        for bond in self.cage.bonds:
            a1id = bond.atom1.id
            a2id = bond.atom2.id
            if a1id in atom_ids_in_G and a2id in atom_ids_in_G:
                cage_g.add_edge(bond.atom1, bond.atom2)

        # Get disconnected subgraphs as molecules.
        # Sort and sort atom ids to ensure molecules are read by RDKIT
        # correctly.
        connected_graphs = [
            sorted(subgraph, key=lambda a: a.id)
            for subgraph in sorted(nx.connected_components(cage_g))
        ]
        for i, cg in enumerate(connected_graphs):
            # Get atoms from nodes.
            atoms = list(cg)
            atom_ids = [i.id for i in atoms]
            sgt = str(len(atoms))
            # Write to mol file.
            filename_ = f'{self.name}_sg{sgt}_{i}.mol'
            self.cage.write(
                filename_,
                atom_ids=atom_ids
            )
            org_lig[filename_] = stk.BuildingBlock.init_from_file(
                filename_
            )
            # Rewrite to fix atom ids.
            org_lig[filename_].write(filename_)
            org_lig[filename_] = stk.BuildingBlock.init_from_file(
                filename_
            )
        return org_lig

    def get_lowest_energy_conformers(self, org_ligs):
        """
        Determine the lowest energy conformer of cage organic linkers.

        Will do multiple if there are multiple types.

        """

        for lig in org_ligs:
            stk_lig = org_ligs[lig]
            # Get optimized ligand name that excludes symmetry.
            opt_lig_n = lig.replace('.mol', '').split('_')[:-1]
            opt_lig_n = opt_lig_n[:-2]+[opt_lig_n[-1]]
            opt_lig_n = '_'.join(opt_lig_n)+'_opt'
            opt_lig_file = f'{opt_lig_n}.mol'

            if not exists(opt_lig_file):
                if not exists(f'{opt_lig_n}_confs/'):
                    os.mkdir(f'{opt_lig_n}_confs/')
                low_e_conf = get_lowest_energy_conformer(
                    name=opt_lig_n,
                    mol=stk_lig
                )
                low_e_conf.write(opt_lig_file)

    def calculate_ligand_SE(self, org_ligs):
        """
        Calculate the strain energy of each ligand in the cage.

        """

        # Check if output file exists.
        if not exists(f'{self.ls_file}.json'):
            strain_energies = {}
            # Iterate over ligands.
            for lig in org_ligs:
                stk_lig = org_ligs[lig]
                ey_file = lig.replace('mol', 'ey')
                # Get optimized ligand name that excludes symmetry.
                opt_lig_n = lig.replace('.mol', '').split('_')[:-1]
                opt_lig_n = opt_lig_n[:-2]+[opt_lig_n[-1]]
                opt_lig_n = '_'.join(opt_lig_n)+'_opt'
                opt_lig_ey = f'{opt_lig_n}.ey'
                opt_lig_file = f'{opt_lig_n}.mol'
                # Calculate energy of extracted ligand.
                if not exists(ey_file):
                    calculate_energy(
                        name=lig.replace('.mol', ''),
                        mol=stk_lig,
                        ey_file=ey_file
                    )
                # Read energy.
                # kJ/mol.
                E_extracted = read_ey(ey_file)

                # Calculate energy of optimised ligand.
                # Load in lowest energy conformer.
                opt_mol = stk.BuildingBlock.init_from_file(
                    opt_lig_file
                )
                if not exists(opt_lig_ey):
                    calculate_energy(
                        name=opt_lig_n,
                        mol=opt_mol,
                        ey_file=opt_lig_ey
                    )
                # Read energy.
                # kJ/mol.
                E_free = read_ey(opt_lig_ey)
                # Add to list the strain energy:
                # (E(extracted) - E(optimised/free))
                lse = E_extracted - E_free
                # kJ/mol.
                strain_energies[lig] = lse

            # Write data.
            with open(f'{self.ls_file}.json', 'w') as f:
                json.dump(strain_energies, f)

        # Get data.
        with open(f'{self.ls_file}.json', 'r') as f:
            strain_energies = json.load(f)
        print('strains', strain_energies)
        return strain_energies

    def calculate_imine_torsions(self, org_ligs):
        """
        Calculate the imine torsion of all ligands in the cage.

        """

        # Iterate over each ligand and find imines.

        # Calculate torsions.

        # Save to list.

        return []

    def calculate_ligand_planarities(self, org_ligs):
        """
        Calculate the change in planarity of the core of all ligands.

        """

        # Iterate over each ligand and find core based on the input
        # molecule and its FGs (i.e. the core is the parts of the
        # input molecule that is not part of the FGs).

        # Calculate planarity of the cores compared to input molecule.

        # Save to list.

        return []

    def analyze_ligand_strain(self, metal_atom_no):
        """
        Analyse cage ligand geometry for strain.

        """

        # Collect the atomic positions of the organic linkers in the
        # cage for analysis.
        org_ligs = self.get_organic_linkers(metal_atom_no)
        self.get_lowest_energy_conformers(org_ligs)

        lse_dict = self.calculate_ligand_SE(org_ligs)
        imine_torsion_list = self.calculate_imine_torsions(org_ligs)
        planarity_list = self.calculate_ligand_planarities(org_ligs)

        self.ls_data = {
            'strain_energies': lse_dict,
            'imine_torsions': imine_torsion_list,
            'core_planarities': planarity_list
        }

    def analyze_metal_strain(self):
        """
        Analyse cage geometry using order parameters.

        """

        # Check if output file exists.
        if not exists(f'{self.op_file}.json'):
            print(f'....analyzing geometry of {self.name}')

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
            os.system('rm temp.xyz')

            # Calculate pore size.
            pw_cage_mol.calculate_pore_diameter_opt()
            pw_cage_mol.calculate_pore_volume_opt()

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


class CageSet:
    """
    Class that builds and analyses stk.ConstructuedMolecules.

    """

    def __init__(
        self,
        name,
        cage_dict,
        complex_dicts,
        ligand_dicts,
        ligand_dir,
        complex_dir
    ):
        self.name = name
        self.properties_file = f'{self.name}_CS.json'
        self.cage_dict = cage_dict
        self.complex_dicts = complex_dicts
        self.ligand_dicts = ligand_dicts
        self.cages_to_build = self.define_cages_to_build(
            ligand_dir,
            complex_dir
        )
        self.built_cage_properties = {}

    def define_cages_to_build(self):
        """
        Defines the name and objects of all cages to build.

        """

        raise NotImplementedError(
            f'Not implemented for {self.__class__}'
        )

    def load_properties(self):
        """
        Load class from JSON file.

        """

        if exists(self.properties_file):
            with open(self.properties_file, 'r') as f:
                self.built_cage_properties = json.load(f)

    def dump_properties(self):
        """
        Dump class to JSON file.

        """
        with open(self.properties_file, 'w') as f:
            json.dump(self.built_cage_properties, f, indent=4)

    def __str__(self):
        return (
            f'{self.__class__.__name__}'
            f'(name={self.name})\n'
            f'{self.cage_dict}'
        )

    def __repr__(self):
        return str(self)

    def _get_no_vertices(self, string):
        """
        Get the number of vertices for a given topology.

        """

        topologies = {
            'm4l4spacer': 8,
            'm8l6face': 14,
            'm6l2l3': 10
        }

        try:
            return topologies[string]
        except KeyError:
            raise KeyError(f'{string} not in {topologies.keys()}')

    def _get_rot_vertices(self, string):
        """
        Get the list of rotatable vertices for a given topology.

        Only ligand vertices are rotatable in this case.

        # TODO: Currently only defined for cube (90 deg). Add tri-face.

        """

        if string in ['m4l4spacer', 'm6l2l3']:
            raise NotImplementedError(
                'Currently only defined for cube (90 deg). Add tri.'
            )

        topologies = {
            'm4l4spacer': [4, 5, 6, 7],
            'm8l6face': [8, 9, 10, 11, 12, 13],
            'm6l2l3': [8, 9, 10]
        }

        try:
            return topologies[string]
        except KeyError:
            raise KeyError(f'{string} not in {topologies.keys()}')

    def _get_ratios(self, n_metals):
        rng = range(0, n_metals+1)
        rats = []
        for i in product(rng, rng):
            if i[0]+i[1] == n_metals:
                rats.append(i)
        return rats

    def _load_complex(self, complex_name, complex_dir):
        complex = stk.BuildingBlock.init_from_file(
            join(complex_dir, f'{complex_name}_opt.mol'),
            functional_groups=[stk.BromoFactory()]
        )

        return complex

    def _load_ligand(self, ligand_name, ligand_dir):
        ligand = stk.BuildingBlock.init_from_file(
            join(ligand_dir, f'{ligand_name}_opt.mol'),
            functional_groups=[stk.BromoFactory()]
        )

        return ligand

    def cage_symmetries(
        self,
        string,
        D_complex,
        L_complex,
        linker,
        check_orientation,
        get_all=False
    ):
        """
        Returns cage symmetries for a given topology.

        """

        no_vertices = self._get_no_vertices(string=string)
        rotatable_vertices = self._get_rot_vertices(string=string)
        # Assumes metal complex vertices is all non-rotatable_vertices.
        complex_vertices = [
            i for i in range(no_vertices)
            if i not in rotatable_vertices
        ]

        if string == 'm4l4spacer':
            symm_list = {}
            raise NotImplementedError(
                'symmetries not defined for m4l4'
            )
        elif string == 'm8l6face':
            symm_list = {}

            if get_all:
                raise NotImplementedError()
                symm_list = symmetries.all_m8l6face_symmetries(
                    D_complex=D_complex,
                    L_complex=L_complex,
                    linker=linker,
                    check_orientation=check_orientation,
                    no_vertices=no_vertices,
                    rotatable_vertices=rotatable_vertices,
                    complex_vertices=complex_vertices
                )

            else:
                # Predefined list of symmetries.
                n_metals = 8
                m8l6_symm = symmetries.M8L6_Symmetry(
                    D_complex=D_complex,
                    L_complex=L_complex,
                    linker=linker,
                    n_metals=n_metals,
                    no_vertices=no_vertices
                )
                symm_list['o1'] = m8l6_symm.o1()
                symm_list['th1'] = m8l6_symm.th1()
                symm_list['th2'] = m8l6_symm.th2()
                symm_list['t1'] = m8l6_symm.t1()
                symm_list['s61'] = m8l6_symm.s61()
                symm_list['s62'] = m8l6_symm.s62()
                symm_list['d31'] = m8l6_symm.d31()
                symm_list['d32'] = m8l6_symm.d32()
                symm_list['c2v'] = m8l6_symm.c2v()
                symm_list['c2h'] = m8l6_symm.c2h()

        elif string == 'm6l2l3':
            symm_list = {}
            raise NotImplementedError(
                'symmetries not defined for m6l2l3'
            )
        else:
            raise KeyError(f'{string} not in defined')

        print(f'{len(symm_list)} symmetries to build')
        return symm_list


class HoCube(CageSet):
    """
    Class that builds and analyses stk.ConstructuedMolecules.

    Represents homoleptic cube cages with all necessary symmetries
    and orientations.

    """

    def __init__(
        self,
        name,
        cage_dict,
        complex_dicts,
        ligand_dicts,
        ligand_dir,
        complex_dir
    ):

        super().__init__(
            name,
            cage_dict,
            complex_dicts,
            ligand_dicts,
            ligand_dir,
            complex_dir
        )

        # Get ligand aspect ratio.
        self.ligand_aspect_ratio = self._get_ligand_AR(ligand_dir)

    def _get_ligand_AR(self, ligand_dir):
        """
        Calculate ligand aspect ratio based on binder positions.

        """

        tet_linker = self._load_ligand(
            ligand_name=self.cage_dict['tetratopic'],
            ligand_dir=ligand_dir
        )
        ligand_AR = calculate_binding_AR(tet_linker)

        return ligand_AR

    def define_cages_to_build(self, ligand_dir, complex_dir):
        """
        Defines the name and objects of all cages to build.

        """

        cages_to_build = []

        # Get Delta and Lambda complexes.
        D_complex_name = [
            i for i in self.complex_dicts if 'del' in i
        ][0]
        L_complex_name = [
            i for i in self.complex_dicts if 'lam' in i
        ][0]
        L_complex = self._load_complex(
            complex_name=L_complex_name,
            complex_dir=complex_dir
        )
        D_complex = self._load_complex(
            complex_name=D_complex_name,
            complex_dir=complex_dir
        )

        D_charge = self.complex_dicts[D_complex_name]['total_charge']
        L_charge = self.complex_dicts[L_complex_name]['total_charge']
        D_free_e = self.complex_dicts[D_complex_name]['unpaired_e']
        L_free_e = self.complex_dicts[L_complex_name]['unpaired_e']

        # Get all linkers and their face orientations.
        tet_prop = self.ligand_dicts[self.cage_dict['tetratopic']]
        tet_linker = self._load_ligand(
            ligand_name=self.cage_dict['tetratopic'],
            ligand_dir=ligand_dir
        )

        # Tetratopic homoleptic cages.
        # Get topology as object to be used in following list.
        tet_topo_name = 'm8l6face'
        tet_topo_fn = available_topologies(string=tet_topo_name)

        symmetries_to_build = self.cage_symmetries(
            string='m8l6face',
            D_complex=D_complex,
            L_complex=L_complex,
            linker=tet_linker,
            check_orientation=tet_prop['check_orientations'],
            get_all=False
        )

        for name_string in symmetries_to_build:
            new_name = (
                f"C_{self.cage_dict['corner_name']}_"
                f"{self.cage_dict['tetratopic']}_"
                f"{name_string}"
            )
            building_blocks = (
                symmetries_to_build[name_string]['building_blocks']
            )
            vertex_alignments = (
                symmetries_to_build[name_string]['vertex_alignments']
            )
            rat = symmetries_to_build[name_string]['ratio']

            # Merge linker and complex charges.
            complex_charge = rat[0]*int(D_charge)
            complex_charge += rat[1]*int(L_charge)
            new_charge = tet_prop['net_charge']*6 + complex_charge

            lig_free_e = [
                int(i)*6 for i in tet_prop['total_unpaired_e']
            ]
            compl_free_e = [
                int(i)*rat[0] + int(j)*rat[1]
                for i, j in zip(D_free_e, L_free_e)
            ]
            new_free_electron_options = []
            for opt in product(lig_free_e, compl_free_e):
                new_free_electron_options.append(opt[0]+opt[1])

            new_cage = Cage(
                name=new_name,
                topology_fn=tet_topo_fn,
                building_blocks=building_blocks,
                vertex_alignments=vertex_alignments,
                topology_string=tet_topo_name,
                charge=new_charge,
                free_electron_options=new_free_electron_options
            )
            cages_to_build.append(new_cage)

        return cages_to_build

    def plot_Y(self, data, ylabel, filename, ylim=None):
        C = '#AFE074'
        M = 'o'

        fig, ax = plt.subplots(figsize=(8, 5))
        x_pos_list = []
        names_list = []
        for i, name in enumerate(data):
            print(i, name)
            X = i+2
            names_list.append(name.split('_')[-1])
            x_pos_list.append(X)
            ax.scatter(
                X,
                data[name],
                c=C,
                edgecolors='k',
                marker=M,
                alpha=1.0,
                s=120
            )
        # Set number of ticks for x-axis
        ax.tick_params(axis='both', which='major', labelsize=16)
        # ax.set_xlabel(r'pore volume [$\mathrm{\AA}^3$]', fontsize=16)
        ax.set_ylabel(ylabel, fontsize=16)
        ax.set_xlim(0, i+3)
        ax.set_ylim(ylim)
        ax.set_xticklabels(names_list)
        ax.set_xticks(x_pos_list)

        fig.tight_layout()
        fig.savefig(
            filename,
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()


class HetPrism(CageSet):
    """
    Class that builds and analyses stk.ConstructuedMolecules.

    Represents heteroleptic prism cages and all necessary homoleptic
    cages.

    """

    def __init__(
        self,
        name,
        cage_dict,
        complex_dicts,
        ligand_dicts,
        ligand_dir,
        complex_dir
    ):

        super().__init__(
            name,
            cage_dict,
            complex_dicts,
            ligand_dicts,
            ligand_dir,
            complex_dir
        )

    def define_cages_to_build(self, ligand_dir, complex_dir):
        """
        Defines the name and objects of all cages to build.

        """

        cages_to_build = []

        # Get Delta and Lambda complex.
        print(self.complex_dicts)
        D_complex_name = [
            i for i in self.complex_dicts if 'del' in i
        ][0]
        L_complex_name = [
            i for i in self.complex_dicts if 'lam' in i
        ][0]
        print(D_complex_name, L_complex_name)
        L_complex = self._load_complex(
            complex_name=L_complex_name,
            complex_dir=complex_dir
        )
        D_complex = self._load_complex(
            complex_name=D_complex_name,
            complex_dir=complex_dir
        )

        D_charge = self.complex_dicts[D_complex_name]['total_charge']
        L_charge = self.complex_dicts[L_complex_name]['total_charge']
        D_free_e = self.complex_dicts[D_complex_name]['unpaired_e']
        L_free_e = self.complex_dicts[L_complex_name]['unpaired_e']
        print(D_charge, D_free_e, L_charge, L_free_e)

        # Get all linkers.
        tet_prop = self.ligand_dicts[self.cage_dict['tetratopic']]
        tet_linker = self._load_ligand(
            ligand_name=self.cage_dict['tetratopic'],
            ligand_dir=ligand_dir
        )
        tri_prop = self.ligand_dicts[self.cage_dict['tritopic']]
        tri_linker = self._load_ligand(
            ligand_name=self.cage_dict['tritopic'],
            ligand_dir=ligand_dir
        )

        print(D_complex, L_complex, tet_linker, tri_linker)

        # Tetratopic homoleptic cages of all symmetries.
        # Get topology as object.
        tet_topo_name, tet_topo = available_topologies(
            string='m8l6face'
        )
        tet_topo = tet_topo()

        tet_n_metals = 8
        tet_ratios = self._get_ratios(tet_n_metals)
        for rat in tet_ratios:
            print(rat)
            new_name = (
                f"C_{self.cage_dict['corner_name']}_"
                f"{self.cage_dict['tetratopic']}_"
                f"{tet_topo_name}_"
                f"d{rat[0]}l{rat[1]}_"
                f"SYMM"
            )
            new_bbs = [D_complex, L_complex, tet_linker]
            new_bb_vertices = {
                D_complex: tet_topo.vertices[:rat[0]],
                L_complex: tet_topo.vertices[rat[0]:rat[0]+rat[1]],
                tet_linker: tet_topo.vertices[tet_n_metals:]
            }
            # Merge linker and complex charges.
            complex_charge = rat[0]*int(D_charge)+rat[1]*int(L_charge)
            new_charge = tet_prop['net_charge']*6 + complex_charge

            print(tet_prop[1])
            lig_free_e = [int(i)*6 for i in tet_prop[1]]
            compl_free_e = [
                int(i)*rat[0] + int(j)*rat[1]
                for i, j in zip(D_free_e, L_free_e)
            ]
            print(lig_free_e, compl_free_e)

            new_free_electron_options = []
            for opt in product(lig_free_e, compl_free_e):
                print(opt)
                new_free_electron_options.append(opt[0]+opt[1])

            print(new_charge, new_free_electron_options)
            new_cage = Cage(
                name=new_name,
                bbs=new_bbs,
                topology=tet_topo,
                topology_string=tet_topo_name,
                bb_vertices=new_bb_vertices,
                charge=new_charge,
                free_electron_options=new_free_electron_options
            )
            cages_to_build.append(new_cage)
            print('NOT BUILDING ALL RATIOS CURRENTLY!!!!!!')
            break

        # Tritopic homoleptic cages of all symmetries.
        # Get topology as object.
        tri_topo_name, tri_topo = available_topologies(
            string='m4l4spacer'
        )
        tri_topo = tri_topo()

        tri_n_metals = 4
        tri_ratios = self._get_ratios(tri_n_metals)
        for rat in tri_ratios:
            print(rat)
            new_name = (
                f"C_{self.cage_dict['corner_name']}_"
                f"{self.cage_dict['tritopic']}_"
                f"{tri_topo_name}_"
                f"d{rat[0]}l{rat[1]}_"
                f"SYMM"
            )
            new_bbs = [D_complex, L_complex, tri_linker]
            new_bb_vertices = {
                D_complex: tri_topo.vertices[:rat[0]],
                L_complex: tri_topo.vertices[rat[0]:rat[0]+rat[1]],
                tri_linker: tri_topo.vertices[tri_n_metals:]
            }
            # Merge linker and complex charges.
            complex_charge = rat[0]*int(D_charge)+rat[1]*int(L_charge)
            new_charge = tri_prop['net_charge']*4 + complex_charge

            print(tri_prop['total_unpaired_e'])
            lig_free_e = [
                int(i)*4 for i in tri_prop['total_unpaired_e']
            ]
            compl_free_e = [
                int(i)*rat[0] + int(j)*rat[1]
                for i, j in zip(D_free_e, L_free_e)
            ]
            print(lig_free_e, compl_free_e)

            new_free_electron_options = []
            for opt in product(lig_free_e, compl_free_e):
                print(opt)
                new_free_electron_options.append(opt[0]+opt[1])

            print(new_charge, new_free_electron_options)
            new_cage = Cage(
                name=new_name,
                bbs=new_bbs,
                topology=tri_topo,
                topology_string=tri_topo_name,
                bb_vertices=new_bb_vertices,
                charge=new_charge,
                free_electron_options=new_free_electron_options
            )
            cages_to_build.append(new_cage)
            print('NOT BUILDING ALL RATIOS CURRENTLY!!!!!!')
            break

        # Prisms of all symmetries.
        # Get topology as object.
        pri_topo_name, pri_topo = available_topologies(
            string='m6l2l3'
        )
        pri_topo = pri_topo()

        pri_n_metals = 6
        pri_ratios = self._get_ratios(pri_n_metals)
        for rat in pri_ratios:
            new_name = (
                f"HeP_{self.cage_dict['corner_name']}_"
                f"{self.cage_dict['tritopic']}_"
                f"{self.cage_dict['tetratopic']}_"
                f"{pri_topo_name}_"
                f"d{rat[0]}l{rat[1]}_"
                f"SYMM"
            )
            new_bbs = [D_complex, L_complex, tri_linker, tet_linker]
            new_bb_vertices = {
                D_complex: pri_topo.vertices[:rat[0]],
                L_complex: pri_topo.vertices[rat[0]:rat[0]+rat[1]],
                tri_linker: pri_topo.vertices[
                    pri_n_metals:pri_n_metals+2
                ],
                tet_linker: pri_topo.vertices[pri_n_metals+2:]
            }
            # Merge linker and complex charges.
            complex_charge = rat[0]*int(D_charge)+rat[1]*int(L_charge)
            new_charge = tri_prop['net_charge']*2
            new_charge += tet_prop['net_charge']*3 + complex_charge

            print(
                tri_prop['total_unpaired_e'],
                tet_prop['total_unpaired_e']
            )
            tri_free_e = [
                int(i)*2 for i in tri_prop['total_unpaired_e']
            ]
            tet_free_e = [
                int(i)*3 for i in tet_prop['total_unpaired_e']
            ]
            compl_free_e = [
                int(i)*rat[0] + int(j)*rat[1]
                for i, j in zip(D_free_e, L_free_e)
            ]
            print(tri_free_e, tet_free_e, compl_free_e)

            new_free_electron_options = []
            for opt in product(tri_free_e, tet_free_e, compl_free_e):
                print(opt)
                new_free_electron_options.append(opt[0]+opt[1]+opt[2])

            print(new_charge, new_free_electron_options)
            new_cage = Cage(
                name=new_name,
                bbs=new_bbs,
                topology=pri_topo,
                topology_string=pri_topo_name,
                bb_vertices=new_bb_vertices,
                charge=new_charge,
                free_electron_options=new_free_electron_options
            )
            cages_to_build.append(new_cage)
            print('NOT BUILDING ALL RATIOS CURRENTLY!!!!!!')
            break

        return cages_to_build

    def plot_min_OPs_avg_PV(self, X, Y, T):
        topo_c_m = {
            'm4l4spacer': ('#E074AF', 'o', r'M$_4$L$_4$'),
            'm8l6face': ('#AFE074', 'X', r'M$_8$L$_6$'),
            'm6l2l3': ('#74AFE0', 'P', r'M$_6$L$^a_2$L$^b_3$')
        }

        Cs = [
            topo_c_m[C.topology_string][0]
            for C in self.cages_to_build
        ]
        print(len(X), len(Y), len(Cs))

        fig, ax = plt.subplots(figsize=(8, 5))
        for x, y, t in zip(X, Y, T):
            ax.scatter(
                x,
                y,
                c=topo_c_m[t][0],
                edgecolors='k',
                marker=topo_c_m[t][1],
                alpha=1.0,
                s=80
            )
        # Set number of ticks for x-axis
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel(r'pore volume [$\mathrm{\AA}^3$]', fontsize=16)
        ax.set_ylabel(r'min. $q_{\mathrm{oct}}$', fontsize=16)
        ax.set_xlim(0, 2000)
        ax.set_ylim(0, 1)

        # Implement legend.
        for i in topo_c_m:
            ax.scatter(
                -1000,
                -1000,
                c=topo_c_m[i][0],
                edgecolors='k',
                marker=topo_c_m[i][1],
                alpha=1.0,
                s=80,
                label=topo_c_m[i][2]
            )

        ax.legend(fontsize=16)

        fig.tight_layout()
        fig.savefig(
            f'{self.name}_minOPsVSporevol.pdf',
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()
