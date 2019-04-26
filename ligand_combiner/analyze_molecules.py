#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module of functions used to analyze all built molecules from DB of core,
ligands and linkers.

Author: Andrew Tarzia

Date Created: 17 Apr 2019

"""

import sys
import matplotlib.pyplot as plt
from matplotlib.cm import RdBu
sys.path.insert(0, '/home/atarzia/thesource/')
from plotting import define_plot_cmap, scatter_plot


def output_analysis(molecule_pop, pair_data, angle_tol, energy_tol):
    '''Output the pair analysis.

    '''
    for p, poly1 in enumerate(molecule_pop):
        combinations = [i for i in pair_data
                        if p == i.popn_ids[0]]
        print(p, poly1.name)
        print(len(combinations))
        if len(combinations) == 0:
            continue
        prefix = poly1.name + '_analysis_'

        print(combinations[0].popn_ids)
        X_data = [max([i.angle1_deviation, i.angle2_deviation])
                  for i in combinations if i.test_N_N_lengths]
        Y_data = [i.NPdN_difference
                  for i in combinations if i.test_N_N_lengths]
        Z_data = [max([i.energy1, i.energy2])/energy_tol
                  for i in combinations if i.test_N_N_lengths]
        # define colour map based on energy tol
        cmap = {'mid_point': energy_tol/2/energy_tol,
                'cmap': RdBu,
                'ticks': [0, energy_tol/2/energy_tol, energy_tol/energy_tol],
                'labels': ['0', str(energy_tol/2), str(energy_tol)],
                'cmap_label': 'energy [kJ/mol]'}

        scatter_plot(X=X_data, Y=Y_data,
                     outfile=prefix+'main.pdf',
                     xtitle='maximum angle deviation [degrees]',
                     ytitle='deviation from planarity',
                     title=poly1.name,
                     xlim=(0, 180), ylim=(0, round(max(Y_data))+1),
                     c='firebrick', edgecolors='k',
                     marker='o', alpha=0.5, s=80, Z=Z_data,
                     cmap=cmap)
        # sys.exit()


def plot_all_pair_info(pair_data, angle_tol, energy_tol):
    '''Do multi dimensional plot of all molecule pair information.

    '''
    markers = ['o', 'X', 'P', 'D', '>', '<', 's', '^', 'p', '*', 'h', 'v']
    fig, ax = plt.subplots(figsize=(8, 5))

    # define colour map based on energy tol
    cmap = define_plot_cmap(fig, ax,
                            mid_point=energy_tol/2/energy_tol, cmap=RdBu,
                            ticks=[0, energy_tol/2/energy_tol, energy_tol/energy_tol],
                            labels=['0', str(energy_tol/2), str(energy_tol)],
                            cmap_label='energy [kJ/mol]')

    X_data = [max([i.angle1_deviation, i.angle2_deviation])
              for i in pair_data if i.test_N_N_lengths]
    Y_data = [i.NPdN_difference
              for i in pair_data if i.test_N_N_lengths]
    Z_data = [max([i.energy1, i.energy2])/energy_tol
              for i in pair_data if i.test_N_N_lengths]
    # M_data = [markers[i.popn_ids[0]]
    #           for i in pair_data if i.test_N_N_lengths]
    X_max = max(X_data)
    Y_max = max(Y_data)

    ax.scatter(X_data, Y_data, c=cmap(Z_data), edgecolors='k',
               marker='o', alpha=0.5, s=80)
    # ax.axhline(pair_data[0].tol, c='k', alpha=0.5)
    # ax.axvline(angle_tol, c='k', alpha=0.5)
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('maximum angle deviation [degrees]', fontsize=16)
    ax.set_ylabel('deviation from planarity', fontsize=16)
    ax.set_xlim(0, round(X_max+1))
    ax.set_ylim(0, round(Y_max+1))
    fig.tight_layout()
    fig.savefig('all_pair_info.pdf', dpi=720,
                bbox_inches='tight')
