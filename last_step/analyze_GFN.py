#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to analyze GFN jobs for all structures.

Author: Andrew Tarzia

Date Created: 08 Mar 2019

"""

import os
import sys
import glob
import json
sys.path.insert(0, '/home/atarzia/thesource/')
import calculations
import GFN_functions


if __name__ == "__main__":
    if (not len(sys.argv) == 2):
        print("""
Usage: analyze_GFN.py dir
    dir: directory to analyze
""")
        sys.exit()
    else:
        targ_dir = sys.argv[1]
    initial_dir = os.getcwd()+'/'
    analysis_dir = initial_dir+targ_dir
    print('analyzing', analysis_dir)
    output_file = targ_dir+'.json'
    initial_struct_dir = initial_dir+'init_structures/'
    xyzs = sorted(glob.glob('*.xyz'))
    # change into dir
    os.chdir(analysis_dir)
    results = {}
    for i in xyzs:
        file = i.replace('.xyz', '')
        # print(file)
        out = file+'.output'
        os.chdir(file+'/')
        # determine properties from GFN output file.
        results[file] = GFN_functions.get_energies(out)
        # calculate RMSD of all structures to input XYZ
        # obviously skip if SPE calculation
        if 'SPE' in targ_dir:
            results[file]['RMSD'] = 0
        else:
            results[file] = calculations.calculate_RMSD(
                                results[file], file,
                                init_structure_dir=initial_struct_dir)
        os.chdir(analysis_dir)
        # print('done')
        # break
    os.chdir(initial_dir)
    # print(results)
    with open(output_file, 'w') as outfile:
        json.dump(results, outfile)
    for i in results:
        print(i, results[i]['TE'])
