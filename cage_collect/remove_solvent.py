#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to remove all non-cage molecules using pyWindow

Author: Andrew Tarzia

Date Created: 15 Apr 2019

"""

import sys
from ase.atoms import Atoms
from os.path import isfile
from ase.visualize import view
sys.path.insert(0, '/home/atarzia/thesource/')
from pywindow_functions import rebuild_system, remove_solvent
from IO_tools import convert_CIF_2_PDB


def main():
    if (not len(sys.argv) == 2):
        print("""
Usage: remove_solvent.py CIF
    CIF (str) - name of CIF to analyze
    """)
        sys.exit()
    else:
        CIF = sys.argv[1]

    if CIF[-4:] != '.cif':
        raise Exception('input file: {} was not a CIF'.format(CIF))
        # no need to redo already done structures
        if isfile(CIF.replace('.cif', '_nosolv.cif')):
            if isfile(CIF.replace('.cif', '_nosolv.pdb')):
                continue

    pdb_file, struct = convert_CIF_2_PDB(CIF)
    # get final struct equivalent to input struct, but without atoms
    final_struct = Atoms()
    final_struct.set_cell(struct.cell)
    final_struct.set_pbc([True, True, True])
    view(struct)
    view(final_struct)
    rebuilt_structure = rebuild_system(file=pdb_file)
    print('modularising...')
    rebuilt_structure.make_modular()
    print('done.')
    # test if one molecule is huge because disorder breaks pywindow code
    no_atoms_orig = len(struct)
    n_atoms_list = []
    for molecule in rebuilt_structure.molecules:
        n_atoms_list.append(rebuilt_structure.molecules[molecule].no_of_atoms)
    max_count = max(n_atoms_list)
    if max_count > no_atoms_orig:
        print('1 UC:', no_atoms_orig, 'modularized max:', max_count)
        # implies that this structure is too disordered for pywindow to handle
        sys.exit('skipping this CIF because modularising failed.')
    final_struct = remove_solvent(pw_struct=rebuilt_structure,
                                  ASE_struct=final_struct,
                                  mol_list=n_atoms_list)
    view(final_struct)
    # output to CIF
    output = CIF.replace('.cif', '_nosolv.cif')
    final_struct.write(output, format='cif')
    output = CIF.replace('.cif', '_nosolv.pdb')
    final_struct.write(output)


if __name__ == "__main__":
    main()
