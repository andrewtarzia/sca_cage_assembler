#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to analyze GFN jobs, where the number of jobs is defined by a built cage
list.

Author: Andrew Tarzia

Date Created: 08 Apr 2019

"""

from os.path import isfile, join
import sys
import json
import glob
sys.path.insert(0, '/home/atarzia/thesource/')
from GFN_functions import get_energies
from stk_functions import topo_2_noimines, expected_window


def write_csv_entry(dict, columns, file):
    line = []
    for c in columns:
        line.append(str(dict[c]))
    with open(file, 'a') as f:
        f.write(','.join(line) + '\n')


def main():
    if (not len(sys.argv) == 3):
        print("""
Usage: analyze_cage_list.py output_file suffix
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
    columns = ['file', 'aldehyde', 'amine', 'topology', 'calc_complete',
               'FE_au', 'no_imines', 'PW_file', 'shape_persist']
    with open(output_file, 'w') as f:
        f.write(','.join(columns) + '\n')
    xyzs = sorted(glob.glob('*' + suffix))
    for i in xyzs:
        cage_prop = {}
        dir = i.rstrip('.xyz')
        gfn_out = i.replace('.xyz', '.output')
        PW_out = i.replace('.xyz', '_properties.json')
        cage_prop['file'] = i
        cage_prop['PW_file'] = PW_out
        # print(i, dir)
        # decompose file name into BBs
        # name convention: aldehyde-name_amine-name_topology
        cage_prop['aldehyde'] = i.split('_')[0]
        cage_prop['amine'] = i.split('_')[1]
        cage_prop['topology'] = i.split('_')[2]
        cage_prop['no_imines'] = topo_2_noimines(cage_prop['topology'])
        # check for GFN output file and collect free energy from GFN output
        # file if available.
        if isfile(join(dir, gfn_out)):
            energies = get_energies(output_file=join(dir, gfn_out))
            try:
                cage_prop['FE_au'] = float(energies['FE'])
                cage_prop['calc_complete'] = 0
            except KeyError:
                cage_prop['FE_au'] = 0
                cage_prop['calc_complete'] = 1
        else:
            cage_prop['FE_au'] = 0
            cage_prop['calc_complete'] = 1
        # get pywindow results
        with open(PW_out, 'r') as f:
            PW_data = json.load(f)
        # set shape persitence flag
        # assume False
        shape_persist = False
        ################################
        # (Computationally-inspired discovery of an unsymmetrical porous
        # organic cage)
        ################################
        # pore diamter >= 3.4 angstrom
        if PW_data['pore_diameter_opt']['diameter'] >= 3.4:
            # structures has >= exp no. windows for topology
            w_no = len(PW_data['windows']['diameters'])
            if w_no >= expected_window(cage_prop['topology']):
                # max window diamter >= 2.8 angstrom
                w_max = max(PW_data['windows']['diameters'])
                if w_max >= 2.8:
                    shape_persist = True
        cage_prop['shape_persist'] = shape_persist
        # output to csv
        # print(cage_prop)
        write_csv_entry(dict=cage_prop, columns=columns, file=output_file)
        # break

    sys.exit()


if __name__ == "__main__":
    main()
