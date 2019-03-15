#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run GFN jobs for all structures.

Author: Andrew Tarzia

Date Created: 08 Mar 2019

"""

import os
import sys
import glob


def setup_dirs(xyzs):
    '''Setup directories for a series of XYZ files (in xyzs).

    '''

    for i in xyzs:
        file = i.replace('.xyz', '')
        print(file)
        os.system('mkdir '+file)
        os.system('cp '+i+' '+file+'/')
        if os.path.isfile('xctrl') is False:
            print('copy an xctrl file in this dir!', os.getcwd())
            sys.exit('exitting.')
        else:
            os.system('cp xctrl '+file+'/')


if __name__ == "__main__":
    xyzs = glob.glob('*.xyz')
    setup_dirs(xyzs)
    # select a GFN execution.
    GFN_exec = '/home/atarzia/software/xtb190301/bin/xtb'
    print('most setup should be in the xctrl file that is in this directory!')
    print('''Option 1: -I xctrl --hess  <<< SPE, no solvent
    Option 2: -I xctrl --hess --gbsa <<< SPE, w solvent in xctrl
    Option 3: -I xctrl --ohess   <<< opt, w solvent in xctrl
    Option 4: -I xctrl --ohess --gbsa  <<< opt, w solvent in xctrl
    ''')
    option = input('select an option!')
    if option == '1':
        part_2 = '-I xctrl --hess >'
    if option == '2':
        part_2 = '-I xctrl --hess --gbsa >'
    if option == '3':
        part_2 = '-I xctrl --ohess >'
    if option == '4':
        part_2 = '-I xctrl --ohess --gbsa >'
    for i in xyzs:
        file = i.replace('.xyz', '')
        print(file)
        out = file+'.output'
        os.chdir(file+'/')
        exec = GFN_exec+' '+i+' '+part_2+' '+out
        os.system(exec)
        os.chdir('../')
        print('done')
