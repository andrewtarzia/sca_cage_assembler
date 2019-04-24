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
    if (not len(sys.argv) == 3):
        print("""
Usage: remove_solvent.py CIF ignore
    CIF (str) - name of CIF to analyze ('*.cif' for all in working dir)
    ignore (str) - string to use to ignore certain files (set NONE if not used)
    """)
        sys.exit()
    else:
        if '*' in sys.argv[1]:
            from glob import glob
            if sys.argv[2] != 'NONE':
                CIFs = sorted([i for i in glob(sys.argv[1]) if sys.argv[2] not in i])
            else:
                CIFs = sorted([i for i in glob(sys.argv[1])])
            print('{} CIFs to analyze'.format(len(CIFs)))
        else:
            CIFs = [sys.argv[1]]

    for CIF in CIFs:
        # no need to redo already done structures
        if isfile(CIF.replace('.cif', '_nosolv.cif')):
            if isfile(CIF.replace('.cif', '_nosolv.pdb')):
                continue
        print('doing', CIF)
        if CIF[-4:] != '.cif':
            raise Exception('input file: {} was not a CIF'.format(CIF))

        pdb_file, struct = convert_CIF_2_PDB(CIF)
        if pdb_file is None and struct is None:
            continue
        # get final struct equivalent to input struct, but without atoms
        final_struct = Atoms()
        final_struct.set_cell(struct.cell)
        final_struct.set_pbc([True, True, True])
        # view(struct)
        # view(final_struct)
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
            # sys.exit('skipping this CIF because modularising failed.')
            print('skipping this CIF because modularising failed.')
            print('----------------------------------------------')
            continue
        final_struct = remove_solvent(pw_struct=rebuilt_structure,
                                      ASE_struct=final_struct,
                                      mol_list=n_atoms_list)
        # only output structures with more than 0 atoms
        if len(final_struct):
            # view(final_struct)
            # output to CIF
            output = CIF.replace('.cif', '_nosolv.cif')
            final_struct.write(output, format='cif')
            output = CIF.replace('.cif', '_nosolv.pdb')
            final_struct.write(output)
        print('done')
        print('----------------------------------------------')


if __name__ == "__main__":
    main()
