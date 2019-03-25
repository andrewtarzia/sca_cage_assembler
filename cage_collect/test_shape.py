#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to get RMSD between two structures.

Author: Andrew Tarzia

Date Created: 08 Mar 2019

"""

import sys
import glob
sys.path.insert(0, '/home/atarzia/thesource/')
import conversion
import pywindow_functions


def main():
    if (not len(sys.argv) == 2):
        print("""
Usage: test_shape.py run_calc
    run_calc: t/T if you want to rerun all pywindow calculations
        (F if just plotting)
    """)
        sys.exit()
    else:
        run_calc = sys.argv[1]

    if run_calc.lower() == 't':
        # get CIFs of interest -- cleaned if it exists, otherwise original
        list_of_cifs = [i for i in glob.glob('*.cif') if '_cleaned' not in i]
        list_of_cleaned_cifs = glob.glob('*_cleaned.cif')
        calculation_list = [i for i in list_of_cifs
                            if i.replace('.cif', '_cleaned.cif') not in list_of_cleaned_cifs]
        for i in list_of_cleaned_cifs:
            calculation_list.append(i)
        print(calculation_list)
        for calc in calculation_list:
            pdb_file, _ = conversion.convert_CIF_2_PDB(calc)
            del _  # we don't need the ASE structure in this case
        # rebuild system
        rebuilt_structure = pywindow_functions.rebuild_system(file=pdb_file)
        rebuilt_structure.make_modular()
        # run analysis on rebuilt system (extracts all cages)
        _ = pywindow_functions.analyze_rebuilt(rebuilt_structure,
                                               file_prefix=pre_op,
                                               atom_limit=20,
                                               include_coms=False,
                                               verbose=False)
        del _  # not needed again


if __name__ == "__main__":
    main()
