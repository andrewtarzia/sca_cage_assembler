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
from numpy.linalg import norm
from numpy import asarray, cos, radians
import matplotlib.pyplot as plt
from matplotlib.cm import RdBu
sys.path.insert(0, '/home/atarzia/thesource/')
from plotting import define_plot_cmap
from Combiner import Combination


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
        print('molecule', poly1.name)
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
            print(poly1.name, poly2.name)
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
                    comb = Combination(poly1, poly2, conf1, conf2)
                    comb.popn_ids = (i, j)
                    # if molecule1 or molecule2 energy > threshold from conf min
                    # skip pair
                    if comb.energy1 > settings['energy_tol'] or comb.energy2 > settings['energy_tol']:
                        # print('skipping conformer pair due to energy')
                        # print(poly1.name, conf1, poly2.name, conf2)
                        # print(comb.energy1, comb.energy2)
                        continue
                    # obtain all properties
                    # check N-N distance of poly1-conf > poly2-conf
                    comb.NN_dist1 = norm(asarray(PROP1['NN_v']))
                    comb.NN_dist2 = norm(asarray(PROP2['NN_v']))
                    # only save combinations with NN_dist1 > NN_dist2
                    # turn on pair specific checks
                    # we dont check for NN distances in this case, i.e.
                    # we DO check if mol_pair is None
                    if mol_pair is None:
                        if comb.test_N_N_lengths() is False:
                            continue
                    # determine angles made by NN_v and NN-BC_v
                    # check that the pairs sum to 180
                    comb.p1_angle1 = PROP1['NN_BCN_1']
                    comb.p1_angle2 = PROP1['NN_BCN_2']
                    comb.p2_angle1 = PROP2['NN_BCN_1']
                    comb.p2_angle2 = PROP2['NN_BCN_2']
                    # now check that the length of the long
                    # vector and the short vector are commensurate
                    # with an ideal trapezoid with the given angles
                    # i.e. the extender vector determined by the
                    # difference of the two NN_v (LHS) matches what is
                    # expected by trig (RHS)
                    comb.extender_V_LHS = (comb.NN_dist1 - comb.NN_dist2) / 2
                    comb.test_angle = radians(180 - comb.p2_angle1)
                    comb.extender_V_RHS = vector_length * cos(comb.test_angle)
                    comb.tol = vector_std * cos(comb.test_angle)
                    # get final geometrical properties
                    comb.get_angle_deviations()
                    comb.get_N_Pd_lengths_deviation()
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
        print(p, poly1.name)
        print(len(combinations))
        if len(combinations) == 0:
            continue
        if mol_pair is None:
            prefix = poly1.name + '_analysis_'
        else:
            prefix = mol_pair[0]+'-'+mol_pair[1]+'_'

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

        fig, ax = plt.subplots(figsize=(8, 5))
        cmp = define_plot_cmap(fig, ax,
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
