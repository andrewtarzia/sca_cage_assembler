#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to evaluate parameter-search DFT calculations.

Author: Andrew Tarzia

Date Created: 20 Jan 2022

"""

import re
import sys
import glob
import numpy as np
import os
import matplotlib.pyplot as plt


def get_cp2k_energy(file):
    search = 'ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:'
    nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
    with open(file, 'r') as f:
        for line in f.readlines():
            if search in line:
                string = nums.search(line.rstrip()).group(0)
                return float(string)


def check_convergence(file):
    nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")

    converge_line = False
    with open(file, 'r') as f:
        for line in f.readlines():
            if '*** SCF run converged in' in line:
                converge_line = True
            if 'Electronic density on regular grids:' in line:
                string = nums.search(line.rstrip()).group(0)
                edense_on_grids = -1*float(string)
            if 'Number of electrons:' in line:
                string = nums.search(line.rstrip()).group(0)
                electrons = float(string)

    if np.isclose(edense_on_grids, electrons, atol=1E-1):
        electron_density = True
    else:
        electron_density = False

    if converge_line and electron_density:
        return True
    else:
        return False


def get_cp2k_forces(file):
    search = 'SUM OF ATOMIC FORCES'
    nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
    with open(file, 'r') as f:
        for line in f.readlines():
            if search in line:
                string = nums.search(line.rstrip()).group(0)
                print(string)
                raise SystemExit()
                return float(string)

    return None


def main():
    if (not len(sys.argv) == 2):
        print(
            """
Usage: evaluate_convergence_tests.py dft_directory

    dft_directory : (str)
        Directory to run dft from - will be created.

""")
        sys.exit()
    else:
        dft_directory = sys.argv[1]

    _figure_path = 'figures'
    output_files = [
        i for i in glob.glob(os.path.join(dft_directory, '*out'))
        if 'slurm' not in i
    ]
    print(len(output_files))

    systems_data = {
        'C_cl1_quad2_12_th2': {},
        'C_cl1_quad2_3_th2': {},
        'C_cl1_quad2_8_s62': {},
    }
    for fi in output_files:
        temp = fi.replace(dft_directory, '').replace('_spe.out', '')
        system = temp.split('_optc_')[0][1:]
        cutoff = int(temp.split('_optc_')[1].split('_')[0])
        rel_cutoff = int(temp.split('_optc_')[1].split('_')[1])
        # Check convergence
        if check_convergence(fi):
            energy = get_cp2k_energy(fi)
            forces = get_cp2k_forces(fi)
            systems_data[system][(cutoff, rel_cutoff)] = (
                energy, forces
            )
        else:
            print(f'{fi} not converged.')
            systems_data[system][(cutoff, rel_cutoff)] = None

    print('cutoff convergence:')
    inital_rel_cutoff = 60
    chosen_cutoff = 700
    for syst in systems_data:
        sdata = systems_data[syst]
        xs = []
        ys = []
        for ds in sdata:
            ey_da = sdata[ds][0]
            fo_da = sdata[ds][1]
            if ey_da is not None and ds[1] == inital_rel_cutoff:
                xs.append(ds[0])
                ys.append(ey_da)

        xys = sorted(zip(xs, ys), key=lambda pair: pair[0])

        fig, ax = plt.subplots(figsize=(8, 5))
        fin_y = [i[1] for i in xys][-1]
        ax.plot(
            [i[0] for i in xys],
            [(i[1]-fin_y)*2625.5 for i in xys],
            lw=3,
            marker='o',
            markersize=12,
            c='firebrick',
        )
        for x, y in xys:
            print(x, y)

        ax.axvline(x=chosen_cutoff, c='k', lw=3)

        # Set number of ticks for x-axis
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel(r'cutoff [Ry]', fontsize=16)
        ax.set_ylabel(r'rel. DFT energy [kJ mol$^{-1}$]', fontsize=16)
        # ax.set_xlim(0, 2500)
        ax.set_ylim(-1, 1)
        fig.tight_layout()
        fig.savefig(
            os.path.join(
                _figure_path, f'{syst}_at{inital_rel_cutoff}rel.pdf'
            ),
            dpi=720,
            bbox_inches='tight',
        )

    raise SystemExit('plot forces too')
    print('\nrel cutoff convergence:')
    chosen_rel_cutoff = 60
    for syst in systems_data:
        sdata = systems_data[syst]
        xs = []
        ys = []
        for ds in sdata:
            ey_da = sdata[ds][0]
            fo_da = sdata[ds][1]
            print(ey_da, fo_da, ds, syst)
            if ey_da is not None and ds[0] == chosen_cutoff:
                xs.append(ds[1])
                ys.append(ey_da)

        xys = sorted(zip(xs, ys), key=lambda pair: pair[0])

        fig, ax = plt.subplots(figsize=(8, 5))
        fin_y = [i[1] for i in xys][-1]
        ax.plot(
            [i[0] for i in xys],
            [(i[1]-fin_y)*2625.5 for i in xys],
            lw=3,
            marker='o',
            markersize=12,
            c='firebrick',
        )
        for x, y in xys:
            print(x, y)


        ax.axvline(x=chosen_rel_cutoff, c='k', lw=3)

        # Set number of ticks for x-axis
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel(r'rel. cutoff []', fontsize=16)
        ax.set_ylabel(r'rel. DFT energy [a.u.]', fontsize=16)
        # ax.set_xlim(0, 2500)
        ax.set_ylim(-1, 1)
        fig.tight_layout()
        fig.savefig(
            os.path.join(
                _figure_path, f'{syst}_at{chosen_rel_cutff}ry.pdf'
            ),
            dpi=720,
            bbox_inches='tight',
        )


if __name__ == "__main__":
    main()
