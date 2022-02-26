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
    output_files = glob.glob(os.path.join(dft_directory, '*out'))
    print(len(output_files))

    cutoffs = range(50, 501, 50)
    rel_cutoffs = range(10, 101, 10)

    systems_data = {
        'C_cl1_quad2_12_t': {},
        'C_cl1_quad2_3_c2v': {},
        'C_cl1_quad2_8_s62': {},
    }
    for fi in output_files:
        temp = fi.replace(dft_directory, '').replace('_spe.out', '')
        system = temp.split('_optc_')[0][1:]
        cutoff = int(temp.split('_optc_')[1].split('_')[0])
        rel_cutoff = int(temp.split('_optc_')[1].split('_')[1])
        energy = get_cp2k_energy(fi)
        systems_data[system][(cutoff, rel_cutoff)] = energy

    for syst in systems_data:
        sdata = systems_data[syst]
        min_energy = min([
            sdata[i] for i in sdata if sdata[i] is not None
        ])
        print(min_energy)
        cu_xs = {i: [] for i in cutoffs}
        cu_ys = {i: [] for i in cutoffs}
        for ds in sdata:
            ey_da = sdata[ds]
            if ey_da is not None:
                # ey_da = (sdata[ds] - min_energy)
                cu_xs[ds[0]].append(ds[1])
                cu_ys[ds[0]].append(ey_da)

        fig, ax = plt.subplots(figsize=(8, 5))
        for cu in cutoffs:
            ax.scatter(
                cu_xs[cu],
                [i-min(cu_ys[cu]) for i in cu_ys[cu]],
                label=f'{cu} Ry',
                s=40,
            )

        # Set number of ticks for x-axis
        ax.set_yscale('log')
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel(r'rel. cutoff []', fontsize=16)
        ax.set_ylabel(r'rel. DFT energy [a.u.]', fontsize=16)
        # ax.set_xlim(0, 2500)
        ax.set_ylim(-0.1, 1)
        ax.legend(fontsize=16)
        fig.tight_layout()
        fig.savefig(
            os.path.join(_figure_path, f'{syst}_convergence.pdf'),
            dpi=720,
            bbox_inches='tight',
        )

    for syst in systems_data:
        sdata = systems_data[syst]
        xs = []
        ys = []
        for ds in sdata:
            ey_da = sdata[ds]
            if ey_da is not None and ds[1] == 60:
                xs.append(ds[0])
                ys.append(ey_da)

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.scatter(
            xs,
            [i-min(ys) for i in ys],
            s=80,
        )
        for x, y in zip(xs, ys):
            print(x, y)

        # Set number of ticks for x-axis
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel(r'cutoff [Ry]', fontsize=16)
        ax.set_ylabel(r'rel. DFT energy [a.u.]', fontsize=16)
        # ax.set_xlim(0, 2500)
        # ax.set_ylim(-0.1, 1)
        fig.tight_layout()
        fig.savefig(
            os.path.join(_figure_path, f'{syst}_at60rel.pdf'),
            dpi=720,
            bbox_inches='tight',
        )


    for syst in systems_data:
        sdata = systems_data[syst]
        xs = []
        ys = []
        for ds in sdata:
            ey_da = sdata[ds]
            if ey_da is not None and ds[0] == 400:
                xs.append(ds[1])
                ys.append(ey_da)

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.scatter(
            xs,
            [i-min(ys) for i in ys],
            s=80,
        )
        for x, y in zip(xs, ys):
            print(x, y)

        # Set number of ticks for x-axis
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel(r'rel. cutoff []', fontsize=16)
        ax.set_ylabel(r'rel. DFT energy [a.u.]', fontsize=16)
        # ax.set_xlim(0, 2500)
        # ax.set_ylim(-0.1, 1)
        fig.tight_layout()
        fig.savefig(
            os.path.join(_figure_path, f'{syst}_at400ry.pdf'),
            dpi=720,
            bbox_inches='tight',
        )


if __name__ == "__main__":
    main()
