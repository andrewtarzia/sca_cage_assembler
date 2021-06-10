#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module defining analysis of crystal structures.

Author: Andrew Tarzia

Date Created: 11 Nov 2020

"""

import numpy as np
import json
import matplotlib.pyplot as plt
from itertools import combinations
import os
import networkx as nx

import stk
import pywindow as pw


from cage import UnexpectedNumLigands
from utilities import (
    convert_symm_names,
    calculate_cube_shape_measure,
    angle_between,
    get_organic_linkers,
    get_atom_distance,
    calculate_abs_imine_torsions,
    calculate_molecule_planarity,
    get_order_values,
    optimize_conformer,
)
import env_set


class XtalCage:
    """
    Generic class that analyses cage structures from x-ray structures.

    """

    def __init__(
        self,
        name,
        pdb_file,
        complex_dicts,
        cage_set_dict
    ):

        self.name = name
        self.pdb_file = pdb_file
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
        os.system('rm temp.xyz')
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

    def collect_lowest_energy_conformer_file(
        self,
        cage_directory,
        n_atoms,
        cage_set,
        ligand_name,
        ligand_directory,
    ):

        # From cage analysis - optimized at solvent level.
        already_run_lowest_energy_cage_filename = (
            f'C_{cage_set}_sg{n_atoms}_1_opt.mol'
        )
        # From flex analysis - not optimized at solvent level.
        already_run_lowest_energy_filename = (
            f'{ligand_directory}/{ligand_name}_loweconf.mol'
        )
        final_filename = f'{self.name}_sg{n_atoms}_1_opt.mol'

        if os.path.exists(already_run_lowest_energy_cage_filename):
            mol = stk.BuildingBlock.init_from_file(os.path.join(
                cage_directory, already_run_lowest_energy_cage_filename
            ))
            mol.write(final_filename)
        else:
            # Load from flex analysis and optimize with solvent.
            mol = stk.BuildingBlock.init_from_file(os.path.join(
                cage_directory, already_run_lowest_energy_filename
            ))
            settings = env_set.crest_conformer_settings(
                solvent=self.cage_set_dict['solvent'],
            ),
            low_e_conf = optimize_conformer(
                name=self.name+'low_e_opt',
                mol=mol,
                opt_level=settings['final_opt_level'],
                charge=settings['charge'],
                no_unpaired_e=settings['no_unpaired_e'],
                max_runs=settings['max_runs'],
                calc_hessian=settings['calc_hessian'],
                solvent=settings['solvent']
            )
            low_e_conf.write(final_filename)

    def get_cage_set_measures(self, cage_directory, cage_set):
        measures_file = os.path.join(
            cage_directory, f'{cage_set}_measures.json'
        )
        with open(measures_file, 'r') as f:
            return json.load(f)

    def plot_Y(self, data, xtal_data, ylabel, filename, ylim=None):
        C = '#AFE074'
        M = 'o'

        fig, ax = plt.subplots(figsize=(8, 5))
        x_pos_list = []
        names_list = []
        for i, name in enumerate(data):
            X = i+2
            names_list.append(convert_symm_names(name.split('_')[-1]))
            x_pos_list.append(X)
            ax.scatter(
                X,
                data[name],
                c=C,
                edgecolors='k',
                marker=M,
                alpha=1.0,
                s=180
            )

        # Add xtal data.
        X = i+1+2
        names_list.append('xtal')
        x_pos_list.append(X)
        ax.scatter(
            X,
            xtal_data,
            c=C,
            edgecolors='k',
            marker='X',
            alpha=1.0,
            s=180
        )
        # Set number of ticks for x-axis
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_ylabel(ylabel, fontsize=16)
        ax.set_xlim(1, i+5)
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

    def get_m_shape(self, mol):

        shapes = calculate_cube_shape_measure(self.name, mol)
        return shapes['CU-8']

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
