#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to setup GFN2-xTB calculations, run them in sequence and collate the
results.

Author: Andrew Tarzia

Date Created: 20 Feb 2019

Directory setup prior to running this script:
---------------------------------------------------------------------------
XYZs:
    - 1(!!) cage with solvent
    - isolated solvent
    - all isolated guests
    - all encapsulated guests
---------------------------------------------------------------------------
Naming conventions:
- contain all XYZ files of structures you want calculated:
    - X_cage_Z_solvent.xyz  <- isolated cage structure (X, Z is an arbitrary)
    - Y_guest_N.xyz <- isolated guest structure of guest N (Y is an arbitrary)
    - cage_guest_N.xyz <- encapsulated guest N
    - Z_solvent.xyz <- isolated solvent

- contain param.txt file with parameters:
    solvent,S     <- solvent to use for continuum
    cage_chrg,C   <- charge on cage
    no_solvent,NS <- number of solvent molecules encapsulated in cage
    no_guest,NG   <- number of guest molecules encapsulated in cage
---------------------------------------------------------------------------
Binding energy equation:
- (20/02/19 - from Felix)
    \deltaG_\mathrm{b} = E_\mathrm{HG} + NS * E_\mathrm{S}
                         - NG * E_\mathrm{G} - E_\mathrm{HS}
    - E_\mathrm{HG} <- free energy of host-guest system
    - E_\mathrm{S}  <- free energy of isolated solvent system
    - E_\mathrm{G}  <- free energy of isolated guest system
    - E_\mathrm{HS} <- free energy of host-solvent system
---------------------------------------------------------------------------

"""

import re
import time
import os
import sys
import glob
from ase.units import Hartree, kJ, mol


def check_directory(cages, solvents, structures, calculations):
    """Check directory conforms.

    """
    # only 1 cage molecule
    if cages > 1 or cages == 0:
        sys.exit('more or less than 1 cage molecules found. Exitting.')
    # only 1 solvent molecule
    if solvents > 1 or solvents == 0:
        sys.exit('more or less than 1 solvent molecules found. Exitting.')
    # check that all guests have there isolated and bound form
    for st in structures:
        print(st)
        print(structures[st]['component'])
        if structures[st]['component'] == 'G':
            # is isolated guest
            GN = str(structures[st]['guest_no'])
            if 'cage_guest_'+GN+'.xyz' not in calculations:
                sys.exit('Missing a host-guest structure. Exitting.')
        elif structures[st]['component'] == 'HG':
            GN = str(structures[st]['guest_no'])
            passed = False
            for calc in calculations:
                print(calc)
                splits = calc.replace('.xyz', '').split("_")
                if calc != structures[st]['file_name'] and len(splits) == 3:
                    if splits[2] == GN:
                        passed = True
            if passed is False:
                sys.exit('Missing an isolated guest structure. Exitting.')


def read_structures():
    """Read in all XYZ structures present and check directory setup.

    """
    structures = {}
    cages = 0
    solvents = 0
    guests = []
    calculations = []
    for file in glob.glob('*.xyz'):
        print(file)
        prefix = file.replace('.xyz', '')
        structures[prefix] = {}
        structures[prefix]['file_name'] = file
        splits = prefix.split("_")
        calculations.append(file)
        if splits[-1] == 'solvent':
            if len(splits) == 2:
                comp = 'S'  # isolated solvent
                name = splits[0]
                guest_no = None
                solvents += 1
            elif len(splits) == 4:
                comp = 'HS'  # cage + solvent
                name = splits[0] + '+' + splits[2]
                guest_no = None
                cages += 1
        elif splits[1] == 'guest':
            if splits[0] == 'cage':
                comp = 'HG'  # cage + guest
                guest_no = int(splits[2])
                name = prefix
            else:
                comp = 'G'  # isolated guest
                name = splits[0]
                guest_no = int(splits[2])
                guests.append(name)
        structures[prefix]['name'] = name
        structures[prefix]['component'] = comp
        structures[prefix]['charged'] = None  # need to determine this later
        structures[prefix]['guest_no'] = guest_no
        print(structures)
        print(structures[prefix])
        # input()
    # check directory contains everything
    print(cages)
    print(solvents)
    print(calculations)
    # input()
    check_directory(cages, solvents, structures, calculations)
    return structures, calculations


def read_param_file(file):
    """Read in parameter file.

    """
    for i in open(file, 'r').readlines():
        print(i)
        item = i.rstrip().split(',')
        print(item)
        if item[0] == 'solvent':
            solv = item[1]
        elif item[0] == 'cage_chrg':
            cage_Q = item[1]
        elif item[0] == 'no_solvent':
            NS = int(item[1])
        elif item[0] == 'no_guest':
            NG = int(item[1])
    return solv, cage_Q, NS, NG


def run_GFN(XYZ, charge=0, solvent=None, hessian=True):
    """Run GFN2-XTB geometry optimization.

    Keyword Arguments:
        XYZ (str) - name of file with XYZ coordinates
        charge (int) - total charge of system
        solvent (str) - solvent to use as continuum. Default is None.
        hessian (bool) - calculate the full hessian (free energies). Default is True.

    """
    # hard code executable
    GFN_exec = 'XTB'
    output = XYZ.replace('.xyz', '.output')
    settings = ''
    if hessian is True:
        settings += '--ohess '
    if charge != 0:
        settings += '--chrg '+str(charge)
    if solvent is not None:
        settings += '--gbsa '+solvent
    execute = GFN_exec+' '+XYZ+' '+settings+' > '+output+' &'
    print(execute)
    os.system(execute)


def check_calculation(filename):
    """Check that the calculation for filename ran successfully.

    Defined by:
        - test1: if 'finished' is in the output file
        - test2: if '*** GEOMETRY OPTIMIZATION CONVERGED' is in the output file
        - test3: if '*** convergence criteria satisfied' is in the output file

    """
    test1 = 'finished'
    pass1 = False
    test2 = '*** GEOMETRY OPTIMIZATION CONVERGED'
    pass2 = False
    test3 = '*** convergence criteria satisfied'
    pass3 = False
    output = filename.replace('.xyz', '.output')
    for line in open(output, 'r'):
        if test1 in line:
            pass1 = True
        if test2 in line:
            pass2 = True
        if test3 in line:
            pass3 = True
    if pass1 is False or pass2 is False or pass3 is False:
        print(filename, 'did not finish.')
        print('pass1:', pass1)
        print('pass2:', pass2)
        print('pass3:', pass3)
        sys.exit('exitting.')
    elif pass1 and pass2 and pass3:
        print('all passed')
        print(pass1, pass2, pass3)


def get_results(results, calc):
    """Get the results of calc.

    Obtained results:
        - free energies <- done
        - existance of imaginary frequencies?

    """
    output = calc+'.output'
    results[calc] = {}
    for line in open(output, 'r'):
        # free energy in a.u.
        if 'TOTAL FREE ENERGY' in line:
            l = line.rstrip().split('<===')
            FE_au = l[0].strip()
            print(FE_au)
            results[calc]['FE_au'] = float(FE_au)
    return results


def convert_au_to_kjmol(value):
    '''Conversion uses ase

    '''
    return value * Hartree * mol / kJ


def analyse_results(structures, results, NS=1, NG=1):
    """Analyse results of calculations.

    Keyword Arguments:

    Analysis:
        - calculate binding energy of all guests

    """
    for st in structures:
        # get constants
        if structures[st]['component'] == 'S':
            # isolated_solvent_energy
            E_S = convert_au_to_kjmol(results[st]['FE_au'])
        if structures[st]['component'] == 'HS':
            # encapsulated solvent energy
            E_HS = convert_au_to_kjmol(results[st]['FE_au'])
    # get properties for each guest
    guests = []
    for st in structures:
        if structures[st]['component'] == 'G' and st not in guests:
            guests.append(st)  # avoid redo of this guest
            # isolated guest energy
            E_G = convert_au_to_kjmol(results[st]['FE_au'])
            guest_no = structures[st]['guest_no']
            # get encapsulated guest energy
            for st_HG in structures:
                if structures[st_HG]['component'] == 'HG':
                    if structures[st_HG]['guest_no'] == guest_no:
                        E_HG = convert_au_to_kjmol(results[st_HG]['FE_au'])
                        break
            print(st_HG, st, guest_no)
            binding_energy =  E_HG + NS * E_S - NG  * E_G - E_HS
            print(st, binding_energy)
            input()
            structures[st]['binding_energy'] = binding_energy

    # output results
    output_results(results, structures)
    # plot all results
    plot_results(results, structures)


def output_results(results, structures):
    """Output the results of all calculations and analysis.

    """

    return None


def plot_results(results, structures):
    """Plot the XX of all calculations and analysis.

    """

    return None



def main():
    """Run script.

    """
    if (not len(sys.argv) == 2):
        print("""
Usage: gfn_binding_energies.py skip_calcs
    skip_calcs: t/T if you want to skip all calculations.""")
        sys.exit()
    else:
        skip = sys.argv[1].lower()
    start_time = time.time()
    wd = os.getcwd()+'/'
    # read in all structures and check for missing items
    structures, calculations = read_structures()
    # read in parameter file
    solv, cage_Q, NS, NG = read_param_file('param.txt')
    print(solv, cage_Q, NS, NG)
    input('param ok?')
    # run all GFN calculations
    if skip != 't':
        # get all necessary properties for to run calculations
        input('up to here!')
        for st in structures:
            print('Running calculation of:', st)
            # check if the calculation is already done

            # make directories for each structure calculation
            DIR = st+'/'
            print(st)
            print(DIR)
            print(wd+DIR)
            input()
            # make this more rigorous and more reliable
            # handle overwriting old XYZ files
            try:
                os.mkdir(wd+DIR)
            except FileExistsError:
                pass
            try:
                os.system('cp '+structures[st]['file_name']+' '+wd+DIR+structures[st]['file_name'])
            except FileExistsError:
                pass
            continue
            os.chdir(wd+DIR)
            print(os.getcwd())
            if structures[st]['charged'] is True:
                pass
            run_GFN(structures[st]['file_name'])
            os.chdir(wd)
            print(os.getcwd())
            # check that the calculation ran successfully
            check_calculation(structures[st]['file_name'])
            print('Finished calculation of:', st)
    # analyse results
    results = {}
    for st in structures:
        print('Running analysis of:', st)
        DIR = st+'/'
        os.chdir(wd+DIR)
        # collect results
        results = get_results(results, st)
        os.chdir(wd)
        print('Finished analysis of:', st)
    analyse_results(structures, results, NS=NS, NG=NG)
    print('--------------------------------------------------------------')
    print('Finished all calculations in', time.time() - start_time, 's')
    print('--------------------------------------------------------------')


if __name__ == "__main__":
    main()
