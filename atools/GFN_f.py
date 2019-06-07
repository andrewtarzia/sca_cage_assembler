#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions that are useful for GFN-xTB usage

Author: Andrew Tarzia

Date Created: 15 Mar 2019
"""
import logging
import os
import re
import sys
import subprocess


def run_GFN():
    '''

    '''


def check_GFN_complete(out_file):
    '''Checks if GFN job is completing by searching out_file for '* finished run on'.

    '''
    lines = open(out_file, 'r').readlines()
    for line in reversed(lines):
        if '* finished run on' in line:
            return True
    return False


def GFN_from_xyzs(xyzs, GFN_exec='/home/atarzia/software/xtb_190418/bin/xtb',
                  P=1):
    '''Run GFN calculation using standard parameters and an xctrl file.

    Parameters
    ----------
    xyzs : :class:`list` of :class:`str`
        List of XYZ files to optimize

    GFN_exec : :class:`str`
        Full path to GFN executable. Defaults to using alias XTB

    P : :class:`int`
        Number of processors to use for calculation. Defaults to 1.

    Returns
    -------
    failed : :class:`list` of :class:`str`
        List of XYZ files whose GFN jobs failed

    '''
    logging.info(f'running GFN using {GFN_exec}')
    if len(xyzs) == 0:
        return []
    setup_dirs(xyzs, xctrl='default')
    # select a GFN execution.
    print('most setup should be in the xctrl file that is in this directory!')
    print('''Option 1: -I xctrl --hess  <<< SPE, no solvent
    Option 2: -I xctrl --hess --gbsa <<< SPE, w solvent in xctrl
    Option 3: -I xctrl --ohess   <<< opt, no solvent in xctrl
    Option 4: -I xctrl --ohess --gbsa  <<< opt, w solvent in xctrl
    ''')
    option = input('select an option!')
    if option == '1':
        part_2 = '-I xctrl --hess '
    if option == '2':
        part_2 = '-I xctrl --hess --gbsa '
    if option == '3':
        part_2 = '-I xctrl --ohess '
    if option == '4':
        part_2 = '-I xctrl --ohess --gbsa '

    # implement number of processors
    part_2 += '-P '+str(P)+' >'

    failed = []
    total = len(xyzs)
    count = 0
    for i in xyzs:
        file = i.replace('.xyz', '')
        logging.info(f'doing {file}, which is {count} out of {total}')
        # move dirs
        os.chdir(file + '/')
        out = file + '.output'
        # check if xtb normal termination has occured
        if check_GFN_complete(out_file=out):
            continue
        exec = GFN_exec + ' ' + i + ' ' + part_2 + ' ' + out
        # run GFN using unlimited stack resources
        returncode = subprocess.call(f'ulimit -s unlimited ; {exec}', shell=True)
        if returncode != 0:
            failed.append(i)
        # return to dir
        os.chdir('../')
        logging.info(f'done')
        count += 1
        sys.exit()

    print('--------------------------------------')
    print('all calculations done. failed XYZs:')
    print(failed)
    return failed


def default_xctrl():
    '''String to write to xctrl file for default GFN calculations.

    '''
    string = '''
$gfn
method=2
$opt
optlevel=2
$scc
temp=300
broydamp=0.4
$thermo
temp=298.15
'''
    return string


def setup_dirs(xyzs, xctrl='default'):
    '''Setup directories for a series of XYZ files (in xyzs).

    '''
    for i in xyzs:
        file = i.replace('.xyz', '')
        print(file)
        if os.path.isdir(file) is False:
            os.system('mkdir ' + file)
        os.system('cp ' + i + ' ' + file + '/')
        if os.path.isfile('xctrl') is False:
            if xctrl == 'default':
                with open('xctrl', 'w') as f:
                    f.write(default_xctrl())
                os.system('cp xctrl ' + file + '/')
            else:
                logging.error(f'copy your xctrl file in this dir!: {os.getcwd()}')
                sys.exit('exitting.')
        else:
            os.system('cp xctrl ' + file + '/')


def get_energies(output_file, GFN_exec):
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
        if GFN_exec == '/home/atarzia/software/xtb_190418/bin/xtb':
            # regex:
            nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
            # free energy in a.u.
            if '| TOTAL FREE ENERGY' in line:
                FE_au = nums.search(line.rstrip()).group(0)
                results_part['FE'] = FE_au
            if '| TOTAL ENERGY' in line:
                TE_au = nums.search(line.rstrip()).group(0)
                results_part['TE'] = TE_au
            # if 'SCC energy    :' in line:
            #     li = line.rstrip().split(':')
            #     SCE_au = li[1].strip()
            #     # print(TE_au)
            #     results_part['SCE'] = SCE_au
            # if 'G(T)' in line:
            #     li = line.rstrip().split(' ')
            #     new_l = [i for i in li if i != '']
            #     # print(new_l)
            #     GT_au = new_l[1].strip()
            #     # print(GT_au)
            #     results_part['GT'] = GT_au
            if '| TOTAL ENTHALPY' in line:
                HT_au = nums.search(line.rstrip()).group(0)
                results_part['HT'] = HT_au
        else:
            # old version
            if 'TOTAL FREE ENERGY' in line:
                # old GFN formatting
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
