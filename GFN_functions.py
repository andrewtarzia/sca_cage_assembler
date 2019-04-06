#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions that are useful for GFN-xTB usage

Author: Andrew Tarzia

Date Created: 15 Mar 2019
"""

import os
import sys


def setup_dirs(xyzs):
    '''Setup directories for a series of XYZ files (in xyzs).

    '''

    for i in xyzs:
        file = i.replace('.xyz', '')
        print(file)
        if os.path.isdir(file) is False:
            os.system('mkdir ' + file)
        os.system('cp ' + i + ' ' + file + '/')
        if os.path.isfile('xctrl') is False:
            print('copy an xctrl file in this dir!', os.getcwd())
            sys.exit('exitting.')
        else:
            os.system('cp xctrl ' + file + '/')


def get_energies(output_file):
    """Get the numbers from output file.

    Obtained results (in a.u.):
        - free energies (FE)
        - absolute energy (TE)
        - SCC energy (SCE)
        - G (GT)
        - H (HT)

    """
    results_part = {}
    for line in open(output_file, 'r'):
        # free energy in a.u.
        if 'TOTAL FREE ENERGY' in line:
            li = line.rstrip().split('<===')
            FE_au = li[0].strip()
            # print(FE_au)
            results_part['FE'] = FE_au
        if 'total E       :' in line:
            li = line.rstrip().split(':')
            TE_au = li[1].strip()
            # print(TE_au)
            results_part['TE'] = TE_au
        if 'SCC energy    :' in line:
            li = line.rstrip().split(':')
            SCE_au = li[1].strip()
            # print(TE_au)
            results_part['SCE'] = SCE_au
        if 'G(T)' in line:
            li = line.rstrip().split(' ')
            new_l = [i for i in li if i != '']
            # print(new_l)
            GT_au = new_l[1].strip()
            # print(GT_au)
            results_part['GT'] = GT_au
        if 'H(T)' in line and 'H(0)-H(T)+PV' not in line:
            li = line.rstrip().split(' ')
            new_l = [i for i in li if i != '']
            # print(new_l)
            HT_au = new_l[1].strip()
            # print(HT_au)
            results_part['HT'] = HT_au
    return results_part
