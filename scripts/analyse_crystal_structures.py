#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to analyse crystal structures from PDB files.

Author: Andrew Tarzia

Date Created: 11 Nov 2020

"""

import sys
from os.path import exists

import stk

from atools import get_organic_linkers

from cage import UnexpectedNumLigands
from cage_set import HoCube
from cage_analysis import analyse_cages, analyse_cage_sets
from utilities import read_lib


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
        # xtal_cage.get_lowest_energy_conformer_file(
        #     cage_directory=cage_directory,
        #     n_atoms=[org_ligs[i].get_num_atoms() for i in org_ligs][0],
        #     cage_set=xtals[xtal]['cage_set'],
        # )

        sys.exit()
        )

        # Full cage analysis.
        cage_data['porediam'] = xtal_cage.get_pore_size()

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

    sys.exit()


if __name__ == "__main__":
    main()
