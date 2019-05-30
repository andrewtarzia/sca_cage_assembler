#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to append window and cage COM's to a CIF.

Author: Andrew Tarzia

Date Created: 19 Feb 2019
"""

import sys
from ase.io import read
import pywindow as pw
import logging
import os
sys.path.insert(0, '/home/atarzia/thesource/')
import pywindow_f


def main():
    if (not len(sys.argv) == 3):
        print("""
        Usage: append_all_COM.py pdb ignore
        pdb: file (.pdb) to analyze and add pseudo atoms to ('*.pdb' for all in working dir)
        ignore (str) - string to use to ignore certain files (set NONE if not used)
        """)
        sys.exit()
    if '*' in sys.argv[1]:
        from glob import glob
        if sys.argv[2] != 'NONE':
            pdbs = sorted([i for i in glob(sys.argv[1]) if sys.argv[2] not in i])
        else:
            pdbs = sorted([i for i in glob(sys.argv[1])])
            logging.info(f'{len(pdbs)} pdbs to analyze')
    else:
        pdbs = [sys.argv[1]]

    count = 1
    for file in pdbs:
        # do not redo
        if os.path.isfile(file.replace('.pdb', '_appended.cif')):
            count += 1
            continue
        logging.info(f'doing {file}: {count} of {len(pdbs)}')
        ASE_structure = read(file)
        if ASE_structure is None:
            count += 1
            continue
        pdb = file
        if '_nosolv' in pdb:
            # if solvent is removed and pdb is used, then this is already the
            # rebuilt structure
            struct = pw.MolecularSystem.load_file(pdb)
            struct.make_modular()
        else:
            # rebuild system
            struct = pywindow_f.modularize(file=pdb)
        # print(struct)
        if struct is None:
            # handle pyWindow failure
            sys.exit(f'pyWindow failure on {pdb}')
        # run analysis
        COM_dict = pywindow_f.analyze_rebuilt(struct,
                                              atom_limit=20,
                                              file_prefix=file.replace('.pdb', ''),
                                              verbose=False, include_coms=True)
        # append atoms to ASE structure as pseudo atoms and write out new CIF
        pywindow_f.append_and_write_COMs(COM_dict, ASE_structure, file, suffix='.pdb')
        count += 1


if __name__ == "__main__":
    # logging.basicConfig(level=logging.DEBUG, format='%(levelname)s-%(message)s')
    # logging.debug(f'Debug mode!')
    logging.basicConfig(level=logging.INFO, format='')
    main()
