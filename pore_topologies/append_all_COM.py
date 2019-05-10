#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to append window and cage COM's to a CIF.

Author: Andrew Tarzia

Date Created: 19 Feb 2019
"""

import sys
from os.path import isfile
sys.path.insert(0, '/home/atarzia/thesource/')
from pywindow_f import append_and_write_COMs, rebuild_system, analyze_rebuilt
from IO_tools import convert_CIF_2_PDB


if __name__ == "__main__":
    if (not len(sys.argv) == 3):
        print("""
Usage: append_all_COM.py CIF ignore
    CIF: file (.cif) to analyze and add pseudo atoms to ('*.cif' for all in working dir)
    ignore (str) - string to use to ignore certain files (set NONE if not used)
    """)
        sys.exit()
    if '*' in sys.argv[1]:
        from glob import glob
        if sys.argv[2] != 'NONE':
            CIFs = sorted([i for i in glob(sys.argv[1]) if sys.argv[2] not in i])
        else:
            CIFs = sorted([i for i in glob(sys.argv[1])])
        print('{} CIFs to analyze'.format(len(CIFs)))
    else:
        CIFs = [sys.argv[1]]

    for file in CIFs:
        # do not redo
        if isfile(file.replace('.cif', '_appended.cif')):
            continue
        pdb_file, ASE_structure = convert_CIF_2_PDB(file)
        if pdb_file is None and ASE_structure is None:
            continue
        # rebuild system
        rebuilt_structure = rebuild_system(file=pdb_file)
        rebuilt_structure.make_modular()
        # run analysis
        COM_dict = analyze_rebuilt(rebuilt_structure, atom_limit=20,
                                   file_prefix=file.replace('.cif', ''),
                                   verbose=False, include_coms=True)
        # append atoms to ASE structure as pseudo atoms and write out new CIF
        append_and_write_COMs(COM_dict, ASE_structure, file)
