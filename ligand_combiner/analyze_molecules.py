#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to analyze all built molecules from DB of core, ligands and linkers.

Author: Andrew Tarzia

Date Created: 17 Apr 2019

"""

import sys
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
