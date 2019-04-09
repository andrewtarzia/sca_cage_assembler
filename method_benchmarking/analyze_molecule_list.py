#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to analyze GFN jobs, where the number of jobs is defined by a list of
molecules.

Author: Andrew Tarzia

Date Created: 09 Apr 2019

"""

from os.path import isfile, join
import sys
import glob
sys.path.insert(0, '/home/atarzia/thesource/')
from GFN_functions import get_energies
from IO_tools import write_csv_entry


def main():
    if (not len(sys.argv) == 3):
        print("""
Usage: analyze_molecule_list.py output_file suffix
    output_file (str) - CSV file to save data to
    suffix (str) - file suffix to analyze GFN calculations of
    """)
        sys.exit()
    else:
        output_file = sys.argv[1]
        suffix = sys.argv[2]
    print('{} will be overwritten!'.format(output_file))
    input('Enter to acknowledge this!')
    # prepare output file
    columns = ['file', 'calc_complete',
               'FE_au', 'TE_au']
    with open(output_file, 'w') as f:
        f.write(','.join(columns) + '\n')
    xyzs = sorted(glob.glob('*' + suffix))
    for i in xyzs:
        prop = {}
        dir = i.rstrip('.xyz')
        gfn_out = i.replace('.xyz', '.output')
        prop['file'] = i
        # check for GFN output file and collect free energy from GFN output
        # file if available.
        if isfile(join(dir, gfn_out)):
            energies = get_energies(output_file=join(dir, gfn_out))
            try:
                prop['FE_au'] = float(energies['FE'])
                prop['TE_au'] = float(energies['TE'])
                prop['calc_complete'] = 0
            except KeyError:
                prop['FE_au'] = 0
                prop['TE_au'] = 0
                prop['calc_complete'] = 1
        else:
            prop['FE_au'] = 0
            prop['TE_au'] = 0
            prop['calc_complete'] = 1
        write_csv_entry(dict=prop, columns=columns, file=output_file)
    sys.exit()


if __name__ == "__main__":
    main()
