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
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
sys.path.insert(0, '/home/atarzia/thesource/')
import Combiner
import plotting


def add_clever_lines(ax):
    '''Add horiz and vertical lines to axes at hardcoded places based on analysis of
    bloch2017 cages.

    '''
    # cage 1
    # max deviation from planarity for all four N-Pd-N vectors
    ax.axhline(y=1.232, c='b', label='cage 1 (GFN)', alpha=0.4, linestyle='-')
    # max angle difference from 180 for all four N-Pd-N vectors
    ax.axvline(x=6.3, c='b', alpha=0.4, linestyle='-')
    # cage 2
    # max deviation from planarity for all four N-Pd-N vectors
    ax.axhline(y=0.701, c='r', label='cage 2 (XRD)', alpha=0.4, linestyle='--')
    # max angle difference from 180 for all four N-Pd-N vectors
    ax.axvline(x=7.1, c='r', alpha=0.4, linestyle='--')


def get_all_pairs(molecule_pop, settings, mol_pair=None):
    '''Get all molecule pairs as ::class::Combination from molecule population.

    '''
    # define bond length vector to use based on N-Pd bond distances extracted
    # from survey
    # N-Pd-N length
    vector_length = 2 * settings['bond_mean']
    vector_std = 2 * settings['bond_std']
    # obtain all pair properties in molecule DB
    # poly1 should be the 'large' molecule, while poly2 should be the 'small'
    # molecule of the pair
    all_pairs = []
    for i, poly1 in enumerate(molecule_pop):
        print(f'molecule: {poly1.name}')
        # turn on pair specific checks
        if mol_pair is not None:
            if poly1.name != mol_pair[0]:
                continue
        for j, poly2 in enumerate(molecule_pop):
            # make sure poly1 != poly2
            if i == j:
                continue
            # turn on pair specific checks
            if mol_pair is not None:
                if poly2.name != mol_pair[1]:
                    continue
            print(f'pair: {poly1.name}, {poly2.name}')
            for conf1 in poly1.geom_prop:
                PROP1 = poly1.geom_prop[conf1]
                # skip conformer if dihedral meant the N's were not on the
                # right side or Ns are pointing the wrong way
                if PROP1['skip'] is True:
                    continue
                for conf2 in poly2.geom_prop:
                    PROP2 = poly2.geom_prop[conf2]
                    # skip conformer if dihedral meant the N's were not on the
                    # right side or Ns are pointing the wrong way
                    if PROP2['skip'] is True:
                        continue
                    # define pair ordering by the NN vector length
                    # mol == larger molecule, smol == smaller molecule.
                    if np.linalg.norm(np.asarray(PROP1['NN_v'])) >= np.linalg.norm(np.asarray(PROP2['NN_v'])):
                        comb = Combiner.Combination(lmol=poly1, smol=poly2,
                                                    lconf=conf1, sconf=conf2)
                        comb.popn_ids = (i, j)
                    else:
                        comb = Combiner.Combination(lmol=poly2, smol=poly1,
                                                    lconf=conf2, sconf=conf1)
                        comb.popn_ids = (j, i)
                    # if molecule1 or molecule2 energy > threshold from conf min
                    # skip pair
                    if comb.lenergy > settings['energy_tol']:
                        continue
                    if comb.senergy > settings['energy_tol']:
                        continue
                    # obtain all properties
                    # check N-N distance of poly1-conf > poly2-conf
                    # only save combinations with lNN_dist > sNN_dist
                    # if mol_pair is not None, turn on pair specific checks
                    # we dont check for NN distances in this case, i.e.
                    # we DO check if mol_pair is None
                    if mol_pair is None:
                        if comb.test_N_N_lengths() is False:
                            continue
                    # get final geometrical properties
                    comb.calculate_planarity(vector_length=vector_length,
                                             vector_std=vector_std)
                    # check that the pairs sum to 180
                    comb.get_angle_deviations()
                    comb.get_planarity_deviation()
                    all_pairs.append(comb)
    return all_pairs


def output_analysis(molecule_pop, pair_data, angle_tol, energy_tol,
                    mol_pair=None):
    '''Output the pair analysis for all pairs with all molecules or just the
    pair in mol_pair.

    '''
    for p, poly1 in enumerate(molecule_pop):
        combinations = [i for i in pair_data
                        if p == i.popn_ids[0]]
        if len(combinations) == 0:
            continue
        if mol_pair is None:
            prefix = poly1.name + '_analysis_'
        else:
            prefix = mol_pair[0]+'-'+mol_pair[1]+'_'
        print(f'molecule: {poly1.name}')
        print(f'no. pairs: {len(combinations)}')
        X_data = [max([i.angle1_deviation, i.angle2_deviation])
                  for i in combinations]
        Y_data = [i.planar_diff
                  for i in combinations]
        Z_data = [max([i.lenergy, i.senergy])/energy_tol
                  for i in combinations]
        # define colour map based on energy tol
        cmap = {'mid_point': energy_tol/2/energy_tol,
                'cmap': cm.RdBu,
                'ticks': [0, energy_tol/2/energy_tol, energy_tol/energy_tol],
                'labels': ['0', str(energy_tol/2), str(energy_tol)],
                'cmap_label': 'energy [kJ/mol]'}

        fig, ax = plt.subplots(figsize=(8, 5))
        cmp = plotting.define_plot_cmap(fig, ax,
                                        mid_point=cmap['mid_point'],
                                        cmap=cmap['cmap'],
                                        ticks=cmap['ticks'],
                                        labels=cmap['labels'],
                                        cmap_label=cmap['cmap_label'])
        ax.scatter(X_data, Y_data, c=cmp(Z_data),
                   edgecolors='k',
                   marker='o', alpha=0.4, s=80)
        # Set number of ticks for x-axis
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('maximum angle deviation [degrees]', fontsize=16)
        ax.set_ylabel('deviation from planarity', fontsize=16)
        ax.set_xlim(0, 180)
        ax.set_ylim(0, round(max(Y_data))+1)
        # add constraint lines
        add_clever_lines(ax)
        if mol_pair is None:
            ax.set_title(poly1.name, fontsize=16)
        else:
            ax.set_title(prefix.replace('_', ' '), fontsize=16)
            ax.legend(fontsize=12)
        fig.tight_layout()
        fig.savefig(prefix+'main.pdf', dpi=720,
                    bbox_inches='tight')
        plt.close()


def plot_all_pair_info(pair_data, angle_tol, energy_tol):
    '''Do multi dimensional plot of all molecule pair information.

    '''
    markers = ['o', 'X', 'P', 'D', '>', '<', 's', '^', 'p', '*', 'h', 'v']
    fig, ax = plt.subplots(figsize=(8, 5))

    # define colour map based on energy tol
    cmap = plotting.define_plot_cmap(fig, ax,
                                     mid_point=energy_tol/2/energy_tol,
                                     cmap=cm.RdBu,
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
