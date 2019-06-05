#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to remove all non-cage molecules using pyWindow

Author: Andrew Tarzia

Date Created: 15 Apr 2019

"""

import sys
import logging
from ase.io import read
from ase.atoms import Atoms
from ase.visualize import view
import os
sys.path.insert(0, '/home/atarzia/thesource/')
import pywindow_f


def main():
    if (not len(sys.argv) == 3):
        print("""
Usage: remove_solvent.py pdb ignore
    pdb (str) - name of pdb to analyze ('*_extracted.pdb' for all in working dir)
    ignore (str) - string to use to ignore certain files (set NONE if not used)
    """)
        sys.exit()
    else:
        if '*' in sys.argv[1]:
            from glob import glob
            if sys.argv[2] != 'NONE':
                pdbs = sorted([i for i in glob(sys.argv[1])
                               if sys.argv[2] not in i and 'nosolv' not in i])
            else:
                pdbs = sorted([i for i in glob(sys.argv[1])])
            print('{} pdbs to analyze'.format(len(pdbs)))
        else:
            pdbs = [sys.argv[1]]

    count = 0
    for pdb in pdbs:
        # no need to redo already done structures
        if os.path.isfile(pdb.replace('.pdb', '_nosolv.cif')):
            if os.path.isfile(pdb.replace('.pdb', '_nosolv.pdb')):
                continue
        logging.info(f'doing {pdb}: {count} of {len(pdbs)}')
        if pdb[-4:] != '.pdb':
            raise Exception(f'input file: {pdb} was not a pdb')

        # pdb_file, struct = IO_tools.convert_CIF_2_PDB(pdb)
        # if pdb_file is None and struct is None:
        #     continue
        struct = read(pdb)
        # get final struct equivalent to input struct, but without atoms
        final_struct = Atoms()
        final_struct.set_cell(struct.cell)
        final_struct.set_pbc([True, True, True])
        # view(struct)
        # view(final_struct)
        rebuilt_structure = pywindow_f.modularize(file=pdb)
        if rebuilt_structure is None:
            # handle pyWindow failure
            sys.exit(f'pyWindow failure on {pdb}')
        # test if one molecule is huge because disorder breaks pywindow code
        no_atoms_orig = len(struct)
        n_atoms_list = []
        for molecule in rebuilt_structure.molecules:
            n_atoms_list.append(rebuilt_structure.molecules[molecule].no_of_atoms)
        max_count = max(n_atoms_list)
        if max_count > no_atoms_orig:
            logging.info(f'1 UC: {no_atoms_orig} modularized max: {max_count}')
            # implies that this structure is too disordered for pywindow to handle
            # sys.exit('skipping this CIF because modularising failed.')
            logging.info(f'skipping this CIF because modularising failed.')
            logging.info(f'----------------------------------------------')
            continue
        final_struct = pywindow_f.remove_solvent(pw_struct=rebuilt_structure,
                                                 ASE_struct=final_struct,
                                                 mol_list=n_atoms_list)
        # only output structures with more than 0 atoms
        if len(final_struct):
            # view(final_struct)
            # output to CIF
            output = pdb.replace('.pdb', '_nosolv.cif')
            final_struct.write(output, format='cif')
            # # turn off PBC and cells for writing pdb
            # final_struct.set_cell([0, 0, 0])
            # final_struct.set_pbc(False)
            output = pdb.replace('.pdb', '_nosolv.pdb')
            final_struct.write(output)
        logging.info(f'done')
        logging.info(f'----------------------------------------------')
        count += 1


if __name__ == "__main__":
    # logging.basicConfig(level=logging.DEBUG, format='%(levelname)s-%(message)s')
    # logging.debug(f'Debug mode!')
    logging.basicConfig(level=logging.INFO, format='')
    main()
