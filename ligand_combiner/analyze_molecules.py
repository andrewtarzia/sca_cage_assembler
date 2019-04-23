#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to analyze all built molecules from DB of core, ligands and linkers.

Author: Andrew Tarzia

Date Created: 17 Apr 2019

"""

import sys
import pickle
from os.path import join
import matplotlib.pyplot as plt


def plot_all_pair_info(pair_data, angle_tol, energy_tol, outfile):
    '''Do multi dimensional plot of all molecule pair information.

    '''
    markers = ['o', 'X', 'P', 'D', '>', '<', ]
    fig, ax = plt.subplots(figsize=(8, 5))

    # define colour map based on energy tol

    X_max = 0
    Y_max = 0
    for pair in pair_data:
        if pair.test_N_N_lengths is False:
            continue
        # print(pair.__dict__)
        pair.get_angle_deviations()
        pair.get_N_Pd_lengths_deviation()
        X = max(pair.angle1_deviation, pair.angle2_deviation)
        Y = pair.NPdN_difference
        Z = max([pair.energy1, pair.energy2])
        X_max = max([X_max, X])
        Y_max = max([Y_max, Y])
        ax.scatter(X, Y, c='firebrick', edgecolors='k',
                   marker=markers[pair.popn_ids[0]], alpha=1.0, s=80)
        # break
    ax.axhline(pair.tol, c='k', alpha=0.5)
    ax.axvline(angle_tol, c='k', alpha=0.5)
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('maximum angle deviation [degrees]', fontsize=16)
    ax.set_ylabel('deviation from planarity', fontsize=16)
    ax.set_xlim(0, round(X_max+1))
    ax.set_ylim(0, round(Y_max+1))
    fig.tight_layout()
    fig.savefig(outfile, dpi=720,
                bbox_inches='tight')


def main():
    if (not len(sys.argv) == 4):
        print("""
    Usage: analyze_molecules.py pair_data angle_tol energy_tol
        pair_data (str) - pickle file with pair data
        angle_tol (float) - tolerance to use on angle matching
        energy_tol (float) - max kJ/mol over min energy conformer to allow
        """)
        sys.exit()
    else:
        pair_data = sys.argv[1]
        angle_tol = float(sys.argv[2])
        energy_tol = float(sys.argv[3])

    proj_dir = '/home/atarzia/projects/ligand_combiner/'
    core_dir = proj_dir + 'cores/'
    liga_dir = proj_dir + 'ligands/'
    link_dir = proj_dir + 'linkers/'
    mole_dir = proj_dir + 'molecules/'

    # load in pair data
    with open(join(mole_dir, pair_data), 'rb') as f:
        all_pairs = pickle.load(f)
    print(len(all_pairs))

    plot_all_pair_info(pair_data=all_pairs,
                       angle_tol=angle_tol, energy_tol=energy_tol,
                       outfile=pair_data.replace('.pkl', '.pdf'))


if __name__ == "__main__":
    main()
