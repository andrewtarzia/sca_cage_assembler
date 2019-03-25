#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions that are useful for converting structure files

Author: Andrew Tarzia

Date Created: 15 Mar 2019
"""

from ase.io import read


def convert_CIF_2_PDB(file):
    '''Convert CIF to PDB file, save and return structure.

    '''
    pdb_file = file.replace('.cif', '.pdb')
    print('converting:', file, 'to', pdb_file)
    structure = read(file)
    # view(structure)
    # input()
    structure.write(pdb_file)
    print('conversion done.')
    return pdb_file, structure


def convert_PDB_2_XYZ(file):
    '''Convert PDB to XYZ file, save and return structure

    '''
    xyz_file = file.replace('.pdb', '.xyz')
    print('converting:', file, 'to', xyz_file)
    structure = read(file)
    # view(structure)
    # input()
    structure.write(xyz_file)
    print('conversion done.')
    return xyz_file, structure
