#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run GFN jobs for all XYZ files in a directory.

Author: Andrew Tarzia

Date Created: 08 Mar 2019

"""

import os
import sys
import glob
sys.path.insert(0, '/home/atarzia/thesource/')
from GFN_functions import setup_dirs


def main():
    if (not len(sys.argv) == 2):
        print("""
Usage: run_GFN.py suffix
    suffix (str) - file suffix to run GFN calculation on
    """)
        sys.exit()
    else:
        suffix = sys.argv[1]
    xyzs = glob.glob('*' + suffix)
    setup_dirs(xyzs)
    # select a GFN execution.
    GFN_exec = '/home/atarzia/software/xtb_190318/bin/xtb'
    print('most setup should be in the xctrl file that is in this directory!')
    print('''Option 1: -I xctrl --hess  <<< SPE, no solvent
    Option 2: -I xctrl --hess --gbsa <<< SPE, w solvent in xctrl
    Option 3: -I xctrl --ohess   <<< opt, w solvent in xctrl
    Option 4: -I xctrl --ohess --gbsa  <<< opt, w solvent in xctrl
    ''')
    option = input('select an option! ')
    if option == '1':
        part_2 = '-I xctrl --hess >'
    if option == '2':
        part_2 = '-I xctrl --hess --gbsa >'
    if option == '3':
        part_2 = '-I xctrl --ohess >'
    if option == '4':
        part_2 = '-I xctrl --ohess --gbsa >'

    failed = []
    total = len(xyzs)
    count = 0
    for i in xyzs:
        file = i.replace('.xyz', '')
        #############
        # add a check for already completed job some how?!
        ############
        print('doing {}, which is {} out of {}'.format(file, count, total))
        out = file + '.output'
        os.chdir(file + '/')
        exec = GFN_exec + ' ' + i + ' ' + part_2 + ' ' + out
        res = os.system(exec)
        if res != 0:
            failed.append(xyzs)
        os.chdir('../')
        print('done')
        count += 1

    print('--------------------------------------')
    print('all calculations done. failed XYZs:')
    print(failed)


if __name__ == "__main__":
    main()
