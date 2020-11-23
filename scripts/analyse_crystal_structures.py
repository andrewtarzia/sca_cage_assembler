#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to analyse crystal structures from PDB files.

Author: Andrew Tarzia

Date Created: 11 Nov 2020

"""

import sys
import numpy as np
import json
from os.path import join
from itertools import combinations
from os import system
import networkx as nx

import stk
import pywindow as pw

from atools import (
    angle_between,
    get_organic_linkers,
    get_atom_distance,
    calculate_ligand_SE,
    calculate_ligand_planarities,
    calculate_abs_imine_torsions,
    calculate_metal_ligand_distance,
    calculate_molecule_planarity,
    get_order_values,
)

from cage import UnexpectedNumLigands
from utilities import read_lib, convert_symm_names, get_plottables


class XtalCage:
    """
    Generic class that analyses cage structures from x-ray structures.

    """

    def __init__(
        self,
        name,
        pdb_file,
        ligand_dict,
        complex_dicts,
        cage_set_dict
    ):

        self.name = name
        self.pdb_file = pdb_file
        self.ligand_dict = ligand_dict
        self.complex_dicts = complex_dicts
        self.cage_set_dict = cage_set_dict

        self.stk_mol = stk.BuildingBlock.init_from_file(pdb_file)
        # Translate to origin.
        self.stk_mol = self.stk_mol.with_centroid([0, 0, 0])
        self.stk_mol.write(f'{name}_stkin.mol')

    def get_metal_atom_nos(self):
        return [i['metal_atom_no'] for i in self.complex_dicts]

    def get_pore_size(self):
        print(f'....analyzing porosity of {self.name}')
        # Load cage into pywindow.
        self.stk_mol.write('temp.xyz')
        pw_cage = pw.MolecularSystem.load_file('temp.xyz')
        pw_cage_mol = pw_cage.system_to_molecule()
        system('rm temp.xyz')
        # Calculate pore size.
        return pw_cage_mol.calculate_pore_diameter_opt()

    def get_organic_linkers(self):
        org_ligs, smiles_keys = get_organic_linkers(
            cage=self.stk_mol,
            metal_atom_nos=(self.get_metal_atom_nos()[0], ),
            file_prefix=f'{self.name}_sg'
        )
        expected_ligands = 1
        num_unique_ligands = len(set(smiles_keys.values()))
        if num_unique_ligands != expected_ligands:
            raise UnexpectedNumLigands(
                f'{self.name} had {num_unique_ligands} unique ligands'
                f', {expected_ligands} were expected. Suggests bad '
                'optimization.'
            )

        return org_ligs, smiles_keys

    def calculate_abs_imine_torsions(self, org_ligs):
        return calculate_abs_imine_torsions(
            org_ligs=org_ligs,
            smarts='[#6]-[#7X2]-[#6X3H1]-[#6X3!H1]',
        )

    def get_lowest_energy_conformer_file(
        self,
        cage_directory,
        n_atoms,
        cage_set
    ):

        already_run_lowest_energy_filename = (
            f'C_{cage_set}_sg{n_atoms}_1_opt.mol'
        )
        print(already_run_lowest_energy_filename)
        mol = stk.BuildingBlock.init_from_file(
            join(cage_directory, already_run_lowest_energy_filename)
        )
        new_filename = f'{self.name}_sg{n_atoms}_1_opt.mol'
        print(new_filename)
        sys.exit()
        mol.write(new_filename)

    def get_cage_set_measures(self, cage_directory, cage_set):
        measures_file = join(
            cage_directory, f'{cage_set}_measures.json'
        )
        with open(measures_file, 'r') as f:
            return json.load(f)
    def get_min_order_value(self):

        op_set = get_order_values(
            mol=self.stk_mol,
            metal=self.get_metal_atom_nos()[0],
            per_site=True
        )
        print(op_set)
        print(op_set.keys())
        target_OPs = [op_set[i]['oct'] for i in op_set]
        print(target_OPs)

        return min(target_OPs)

    def define_faces(self, m_structure):

        def update_connections(connections, id1, id2, distance):

            if len(connections[id1]) < 3:
                connections[id1].append((id2, distance))
                connections[id1].sort(key=lambda x: x[1])
            else:
                old_list = connections[id1]
                if distance < max([i[1] for i in old_list]):
                    old_list.pop(-1)
                    new_list = old_list+[(id2, distance)]
                    new_list.sort(key=lambda x: x[1])
                else:
                    new_list = old_list.copy()

                connections[id1] = new_list

            return connections

        def get_connections(structure):

            connections = {
                i.get_id(): [] for i in structure.get_atoms()
            }
            for m_pair in combinations(structure.get_atoms(), 2):
                m1_id = m_pair[0].get_id()
                m2_id = m_pair[1].get_id()
                distance = get_atom_distance(structure, m1_id, m2_id)

                # Add to connections.
                # If length of connections is < 3, add.
                # Else, only update and sort if distance is shorter
                # than an existing one.
                update_connections(connections, m1_id, m2_id, distance)
                update_connections(connections, m2_id, m1_id, distance)

            return connections

        # Find distinct faces by determining which metals are close.
        m_connections = get_connections(m_structure)

        # Use graph cycles to find faces.
        atoms = [i for i in m_structure.get_atoms()]
        metal_graph = nx.Graph()
        in_graph = set()
        for metal in m_connections:
            if metal not in in_graph:
                metal_graph.add_node(atoms[metal])
                in_graph.add(metal)
            for m2, distance in m_connections[metal]:
                if m2 not in in_graph:
                    metal_graph.add_node(atoms[m2])
                    in_graph.add(m2)
                metal_graph.add_edge(atoms[metal], atoms[m2])

        cycles = nx.simple_cycles(metal_graph.to_directed())
        face_atoms = [i for i in cycles if len(i) == 4]
        filtered_face_atoms = []
        face_atom_set = set()
        for fa in face_atoms:
            ordered_fa = tuple(sorted([i.get_id() for i in fa]))
            if ordered_fa not in face_atom_set:
                filtered_face_atoms.append(fa)
                face_atom_set.add(ordered_fa)
        if len(filtered_face_atoms) != 6:
            raise ValueError(
                f'{len(filtered_face_atoms)} faces found, expected 6!'
            )

        self.faces = {}
        for i, fa in enumerate(filtered_face_atoms):
            for j, fa2 in enumerate(filtered_face_atoms):
                not_in_fa = [i for i in fa2 if i not in fa]
                if len(not_in_fa) == 4:
                    opposite_id = j
            self.faces[i] = (fa, opposite_id)

    def get_max_face_metal_PD(self, mol):

        plane_devs = []
        for face in self.faces:
            atom_ids = [i.get_id() for i in self.faces[face][0]]
            plane_devs.append(calculate_molecule_planarity(
                mol=mol,
                plane_ids=atom_ids,
                atom_ids=atom_ids,
            ))
        max_face_metal_PD = max(plane_devs)
        return max_face_metal_PD

    def get_max_face_interior_angle_dev(self, mol):

        pos_mat = mol.get_position_matrix()
        sum_interior_angles = []
        for face in self.faces:
            interior_angles = []
            atom_ids = [i.get_id() for i in self.faces[face][0]]*2
            angles = [
                atom_ids[i: i + 3]
                for i in range(0, len(atom_ids))
            ][: 4]
            for trio in angles:
                vector1 = pos_mat[trio[1]]-pos_mat[trio[0]]
                vector2 = pos_mat[trio[1]]-pos_mat[trio[2]]
                interior_angles.append(np.degrees(
                    angle_between(vector1, vector2)
                ))
            sum_interior_angles.append(sum(interior_angles))
        max_face_interior_angle_dev = max(
            [abs(360-i) for i in sum_interior_angles]
        )
        return max_face_interior_angle_dev

    def get_max_face_anisotropy(self, mol):

        face_anisos = {}
        for face in self.faces:
            atom_ids = [i.get_id() for i in self.faces[face][0]]

            # Pick one metal - use index 1 as central atom.
            d1 = get_atom_distance(
                molecule=mol,
                atom1_id=atom_ids[0],
                atom2_id=atom_ids[1],
            )
            d2 = get_atom_distance(
                molecule=mol,
                atom1_id=atom_ids[1],
                atom2_id=atom_ids[2],
            )
            face_anisos[face] = d2 / d1 if d2 > d1 else d1 / d2
        paired_face_anisos = [
            (i, j, face_anisos[i], face_anisos[j])
            for i, j in combinations(face_anisos, r=2)
            if self.faces[i][1] == j
        ]
        max_face_aniso_diff = max([
            100*((i[2] - i[3]) / i[2]) for i in paired_face_anisos
        ])
        return max_face_aniso_diff

    def write_metal_atom_structure(self):

        metal_atom_ids = [
            i.get_id() for i in self.stk_mol.get_atoms()
            if i.get_atomic_number() in self.get_metal_atom_nos()
        ]

        # Write to mol file.
        self.stk_mol.write(
            f'{self.name}_M.mol',
            atom_ids=metal_atom_ids,
        )

    def __str__(self):
        return (
            f'{self.__class__.__name__}'
            f'(name={self.name})'
        )

    def __repr__(self):
        return str(self)


def main():
    first_line = (
        'Usage: analyse_crystal_structures.py lig_lib_file '
        'prism_lib_file compl_lib_file lig_directory compl_directory '
        'cage_directory'
    )
    if (not len(sys.argv) == 7):
        print(f"""
{first_line}

    ligand_lib_file : (str)
        File containing ligand information (XXXXX)

    complex_lib_file : (str)
        File containing complex information (XXXXX).

    cage_set_lib_file : (str)
        File containing cage information (XXXXX).

    ligand_directory : (str)
        Directory with required ligand structures.

    complex_directory : (str)
        Directory with required complex structures.

    cage_directory : (str)
        Directory with required cage structures.

    """)
        sys.exit()
    else:
        ligand_lib_file = sys.argv[1]
        complex_lib_file = sys.argv[2]
        cage_set_lib_file = sys.argv[3]
        ligand_directory = sys.argv[4]
        complex_directory = sys.argv[5]
        cage_directory = sys.argv[6]

    cage_set_lib = read_lib(cage_set_lib_file)
    complexes = read_lib(complex_lib_file)
    ligands = read_lib(ligand_lib_file)

    # List of the xtal structures and their corresponding names.
    xtals = {
        'jd235': {
            'cage_set': 'cl1_quad2_12',
            'cage_name': 'C_cl1_quad2_12_th1',  # or th2
            'symmetry_name': 'th1',  # or th2
            'ligand_name': 'quad2_12',
            'complexes': ('cl1_zn_oct_lam', 'cl1_zn_oct_del'),
        },
        'jd257': {
            'cage_set': 'cl1_quad2_8',
            'cage_name': 'C_cl1_quad2_8_th1',  # or th2
            'symmetry_name': 'th1',  # or th2
            'ligand_name': 'quad2_8',
            'complexes': ('cl1_zn_oct_lam', 'cl1_zn_oct_del'),
        },
        'jd301': {
            'cage_set': 'cl1_quad2_3',
            'cage_name': 'C_cl1_quad2_3_th1',  # or th2
            'symmetry_name': 'th1',  # or th2
            'ligand_name': 'quad2_3',
            'complexes': ('cl1_zn_oct_lam', 'cl1_zn_oct_del'),
        },
        'jd326': {
            'cage_set': 'cl1_quad2_16',
            'cage_name': 'C_cl1_quad2_16_s61',  # or s62
            'symmetry_name': 's61',  # or s62
            'ligand_name': 'quad2_16',
            'complexes': ('cl1_zn_oct_lam', 'cl1_zn_oct_del'),
        },
    }

    xtal_cage_data = {}
    for xtal in xtals:
        print(f'---- doing: {xtal}')
        pdb_file = f'{xtal}.pdb'
        cage_data = {}
        xtal_cage = XtalCage(
            name=xtal,
            pdb_file=pdb_file,
            ligand_dict=ligands[xtals[xtal]['ligand_name']],
            complex_dicts=[
                complexes[i] for i in xtals[xtal]['complexes']
            ],
            cage_set_dict=cage_set_lib[xtals[xtal]['cage_set']]
        )
        print(xtal_cage)
        org_ligs, smiles_keys = xtal_cage.get_organic_linkers()
        xtal_cage.write_metal_atom_structure()
        xtal_cage.get_lowest_energy_conformer_file(
            cage_directory=cage_directory,
            n_atoms=[org_ligs[i].get_num_atoms() for i in org_ligs][0],
            cage_set=xtals[xtal]['cage_set'],
        )

        # Face-based analysis.
        m_structure = stk.BuildingBlock.init_from_file(
            f'{xtal_cage.name}_M.mol'
        )
        xtal_cage.define_faces(m_structure)
        cage_data['maxfacemetalpd'] = xtal_cage.get_max_face_metal_PD(
            m_structure
        )
        cage_data['maxintangledev'] = (
            xtal_cage.get_max_face_interior_angle_dev(m_structure)
        )
        cage_data['maxdifffaceaniso'] = (
            xtal_cage.get_max_face_anisotropy(m_structure)
        )

        # Full cage analysis.
        cage_data['porediam'] = xtal_cage.get_pore_size()
        cage_data['ML_lengths'] = calculate_metal_ligand_distance(
            mol=xtal_cage.stk_mol,
            metal_atomic_number=30,
            ligand_atomic_number=7,
        )
        cage_data['maxMLlength'] = max(cage_data['ML_lengths'])
        cage_data['octop'] = xtal_cage.get_min_order_value()

        # Ligand analysis.
        cage_data['core_planarities'] = calculate_ligand_planarities(
            org_ligs=org_ligs
        )
        cage_data['imine_torsions'] = (
            xtal_cage.calculate_abs_imine_torsions(org_ligs)
        )
        cage_data['strain_energies'] = calculate_ligand_SE(
            org_ligs=org_ligs,
            smiles_keys=smiles_keys,
            output_json=f'{xtal_cage.name}_lse.json',
            file_prefix=f'{xtal_cage.name}_sg'
        )
        cage_data['lsesum'] = sum([
            cage_data['strain_energies'][i]
            for i in cage_data['strain_energies']
        ])
        cage_data['minitors'] = min([
            j for i in cage_data['imine_torsions']
            for j in cage_data['imine_torsions'][i]
        ])
        cage_data['maxcrplan'] = max([
            cage_data['core_planarities'][i]
            for i in cage_data['core_planarities']
        ])
        xtal_cage_data[xtal] = cage_data
        comp_cage_data[xtal] = xtal_cage.get_cage_set_measures(
            cage_directory=cage_directory,
            cage_set=xtals[xtal]['cage_set'],
        )

        plottables = get_plottables(
            measures=comp_cage_data[xtal],
            name=xtal_cage.name
        )


    sys.exit()


if __name__ == "__main__":
    main()
