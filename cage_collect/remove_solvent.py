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

    pdb_file, struct = convert_CIF_2_PDB(CIF)
    # get final struct equivalent to input struct, but without atoms
    final_struct = Atoms()
    final_struct.set_cell(struct.cell)
    final_struct.set_pbc([True, True, True])
    view(struct)
    view(final_struct)
    rebuilt_structure = rebuild_system(file=pdb_file)
    rebuilt_structure.make_modular()
    final_struct = remove_solvent(pw_struct=rebuilt_structure,
                                  ASE_struct=final_struct)


if __name__ == "__main__":
    main()
