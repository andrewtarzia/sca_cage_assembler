#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to append window and cage COM's to a CIF.

Author: Andrew Tarzia

Date Created: 19 Feb 2019
"""

import sys
import ase
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

    for file in pdbs:
        # do not redo
        if os.path.isfile(file.replace('.pdb', '_appended.cif')):
            continue
        ASE_structure = ase.io.read(file)
        if ASE_structure is None:
            continue
        # rebuild system
        pdb = file
        rebuilt_structure = pywindow_f.modularize(file=pdb)
        if rebuilt_structure is None:
            # handle pyWindow failure
            sys.exit(f'pyWindow failure on {pdb}')
        # run analysis
        COM_dict = pywindow_f.analyze_rebuilt(rebuilt_structure,
                                              atom_limit=20,
                                              file_prefix=file.replace('.pdb', ''),
                                              verbose=False, include_coms=True)
        # append atoms to ASE structure as pseudo atoms and write out new CIF
        pywindow_f.append_and_write_COMs(COM_dict, ASE_structure, file, suffix='.pdb')


if __name__ == "__main__":
    # logging.basicConfig(level=logging.DEBUG, format='%(levelname)s-%(message)s')
    # logging.debug(f'Debug mode!')
    logging.basicConfig(level=logging.INFO, format='')
    main()
