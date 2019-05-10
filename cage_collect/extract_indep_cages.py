#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to extract all independant cages from a CIF using pyWindow code.

Author: Andrew Tarzia

Date Created: 04 Apr 2019

"""

import sys
sys.path.insert(0, '/home/atarzia/thesource/')
import pywindow_f
import IO_tools


def main():
    if (not len(sys.argv) == 2):
        print("""
Usage: extract_indep_cages.py CIF
    CIF (str) - name of CIF to analyze
    """)
        sys.exit()
    else:
        CIF = sys.argv[1]

    if CIF[-4:] != '.cif':
        raise Exception('input file: {} was not a CIF'.format(CIF))

    pdb_file, struct = IO_tools.convert_CIF_2_PDB(CIF)
    if pdb_file is None and struct is None:
        sys.exit()
    rebuilt_structure = pywindow_f.rebuild_system(file=pdb_file)
    rebuilt_structure.make_modular()
    res = pywindow_f.analyze_rebuilt(rebuilt_structure, file_prefix=CIF.rstrip('.cif'),
                                     atom_limit=20, include_coms=False, verbose=False)
    print('===================================================')
    print('Results of pyWindow analysis on all indep cages:')
    print('===================================================')
    for i in res:
        try:
            print('cage {}:'.format(i))
            print('has {} windows of diameter:'.format(len(res[i][1]['diameters'])))
            print(res[i][1]['diameters'])
        except TypeError:
            print(res[i])


if __name__ == "__main__":
    main()
