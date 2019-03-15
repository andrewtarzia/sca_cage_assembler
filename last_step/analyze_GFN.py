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
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd


def calculate_RMSD(results_part, structure_name, init_structure_dir):
    '''Calculate RMSD of GFN output structure compared to initial structure.

    Code from James Pegg.

    '''
    # read in initial structure
    initial_structure_file = init_structure_dir+file+'.xyz'
    ref = mda.Universe(initial_structure_file)
    # read in new structure
    new_structure_file = 'xtbopt.xyz'
    mobile = mda.Universe(new_structure_file)
    # RMSD of all atoms
    RMSD = rmsd(mobile.atoms.positions,
                ref.atoms.positions,
                center=True, superposition=True)
    results_part['RMSD'] = RMSD
    return results_part


def get_results(results_part, output_file):
    """Get the numbers from output file.

    Obtained results (in a.u.):
        - free energies (FE)
        - absolute energy (TE)
        - SCC energy (SCE)
        - deltaG (GT)
        - deltaH (HT)

    """
    for line in open(output_file, 'r'):
        # free energy in a.u.
        if 'TOTAL FREE ENERGY' in line:
            l = line.rstrip().split('<===')
            FE_au = l[0].strip()
            # print(FE_au)
            results_part['FE'] = FE_au
        if 'total E       :' in line:
            l = line.rstrip().split(':')
            TE_au = l[1].strip()
            # print(TE_au)
            results_part['TE'] = TE_au
        if 'SCC energy    :' in line:
            l = line.rstrip().split(':')
            SCE_au = l[1].strip()
            # print(TE_au)
            results_part['SCE'] = SCE_au
        if 'G(T)' in line:
            l = line.rstrip().split(' ')
            new_l = [i for i in l if i != '']
            # print(new_l)
            GT_au = new_l[1].strip()
            # print(GT_au)
            results_part['GT'] = GT_au
        if 'H(T)' in line and 'H(0)-H(T)+PV' not in line:
            l = line.rstrip().split(' ')
            new_l = [i for i in l if i != '']
            # print(new_l)
            HT_au = new_l[1].strip()
            # print(HT_au)
            results_part['HT'] = HT_au
    return results_part


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
        results[file] = {}
        out = file+'.output'
        os.chdir(file+'/')
        # determine properties from GFN output file.
        results[file] = get_results(results[file], out)
        # calculate RMSD of all structures to input XYZ
        # obviously skip if SPE calculation
        if 'SPE' in targ_dir:
            results[file]['RMSD'] = 0
        else:
            results[file] = calculate_RMSD(results[file], file,
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
