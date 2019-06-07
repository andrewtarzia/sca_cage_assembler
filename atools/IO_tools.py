#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions that are useful for reading/writing structure files

Author: Andrew Tarzia

Date Created: 15 Mar 2019
"""
from os.path import isfile
from stk import StructUnit, OPTIONS
from ase.io import read
from ase.io.xyz import write_xyz
from pymatgen.io.cif import CifParser


def convert_MOL3000_2_PDB_XYZ(file):
    '''Convert MOL from stk to PDB and XYZ file. Return None.

    '''
    OPTIONS['cache'] = False  # turn caching off for loading
    pdbfile = file.replace('.mol', '.pdb')
    # read in mol into stk
    struct = StructUnit(file)
    struct.write(pdbfile)
    # read pdb into ASE
    convert_PDB_2_XYZ(file=pdbfile)


def check_ASE_handle(file, wstruct=True):
    '''Check if ASE handles the reading of a given file.

    '''
    try:
        structure = read(file)
        if wstruct:
            return file, structure
        else:
            return file
    except IndexError:
        print('ASE load failed with IndexError. Skipping...')
        if wstruct:
            return None, None
        else:
            return None
    except ValueError:
        print('ASE load failed with IndexError. Skipping...')
        if wstruct:
            return None, None
        else:
            return None


def convert_CIF_2_PDB(file, wstruct=True):
    '''Convert CIF to PDB file, save and return structure.

    '''
    pdb_file = file.replace('.cif', '.pdb')
    print('converting:', file, 'to', pdb_file)
    if isfile(pdb_file) is False:
        try:
            structure = read(file)
        except IndexError:
            print('ASE load failed with IndexError. Skipping...')
            if wstruct:
                return None, None
            else:
                return None
        except ValueError:
            print('ASE load failed with IndexError. Skipping...')
            if wstruct:
                return None, None
            else:
                return None
        structure.write(pdb_file)
        print('conversion done.')
    structure = read(pdb_file)
    if wstruct:
        return pdb_file, structure
    else:
        return pdb_file


def convert_PDB_2_XYZ(file, comment=None):
    '''Convert PDB to standard (NOT exteneded) XYZ file, save and
    return structure

    '''
    xyz_file = file.replace('.pdb', '.xyz')
    print('converting:', file, 'to', xyz_file)
    if isfile(xyz_file) is False:
        structure = read(file)
        if comment is None:
            cmt = 'This is an XYZ structure.'
        else:
            cmt = comment
        write_xyz(xyz_file, images=structure, comment=cmt,
                  columns=['symbols', 'positions'])
        print('conversion done.')
    structure = read(xyz_file)
    return xyz_file, structure


def read_cif_pmg(file, primitive=False):
    '''A function to read CIFs with pymatgen and suppress warnings.

    '''
    s = CifParser(file, occupancy_tolerance=100)
    struct = s.get_structures(primitive=primitive)[0]
    return struct


def write_csv_entry(dict, columns, file):
    line = []
    for c in columns:
        line.append(str(dict[c]))
    with open(file, 'a') as f:
        f.write(','.join(line) + '\n')
