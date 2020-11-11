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


    def __str__(self):
        return (
            f'{self.__class__.__name__}'
            f'(name={self.name})'
        )

    def __repr__(self):
        return str(self)


def get_metal_atom_structure(cage, metal_atom_nos, file_prefix=None):
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

    metal_atom_ids = [
        i.get_id() for i in cage.get_atoms()
        if i.get_atomic_number() in metal_atom_nos
    ]
    print(metal_atom_ids)

    # Write to mol file.
    cage.write(f'{file_prefix}.mol', atom_ids=metal_atom_ids)


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

    for xtal in xtals:
        print(f'---- doing: {xtal}')
        pdb_file = f'{xtal}.pdb'
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
        xtal_cage.get_metal_atom_no
        xtal_cage.get_organic_linkers
        xtal_cage.get_pore_size
        xtal_cage.get_metal_atom_structure
        sys.exit()
        # Ligand strain.
        metal_atom_no = [i['metal_atom_no'] for i in complex_dicts][0]
        org_ligs, smiles_keys = get_organic_linkers(
            cage=stk_mol,
            metal_atom_nos=(metal_atom_no, ),
            file_prefix=f'{xtal}_sg'
        )
        expected_ligands = 6
        num_unique_ligands = len(set(smiles_keys.values()))
        if num_unique_ligands != expected_ligands:
            raise UnexpectedNumLigands(
                f'{xtal} had {num_unique_ligands} unique ligands'
                f', {expected_ligands} were expected. Suggests bad '
                'optimization.'
            )

        # Porosity analysis.

        # Metal-based analysis.
        get_metal_atom_structure(
            cage=stk_mol,
            metal_atom_nos=(metal_atom_no, ),
            file_prefix=f'{xtal}_M'
        )

    sys.exit()


if __name__ == "__main__":
    main()
