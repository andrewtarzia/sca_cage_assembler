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
import numpy as np
import pandas as pd
from rdkit.Chem import AllChem as Chem
from os.path import join
import matplotlib.pyplot as plt
import stk
from stk import Population
from Combiner import get_molecule, get_geometrical_properties, Combination
sys.path.insert(0, '/home/atarzia/thesource/')
from stk_functions import build_population
from rdkit_functions import mol_list2grid


def main():
    if (not len(sys.argv) == 3):
        print("""
    Usage: analyze_molecules.py pair_data angle_tol
        pair_data (str) - pickle file with pair data
        angle_tol (float) - tolerance to use on angle matching
        """)
        sys.exit()
    else:
        pair_data = sys.argv[1]
        angle_tol = sys.argv[2]

    proj_dir = '/home/atarzia/projects/ligand_combiner/'
    core_dir = proj_dir + 'cores/'
    liga_dir = proj_dir + 'ligands/'
    link_dir = proj_dir + 'linkers/'
    mole_dir = proj_dir + 'molecules/'

    # load in pair data
    with open(join(mole_dir, pair_data), 'rb') as f:
        all_pairs = pickle.load(f)
    print(len(all_pairs))
                #
                #
                #
                # NN_dist2s = []
                # angles1_sums = []
                # angles2_sums = []
                # LHSs = []
                # RHSs = []
                #
                #
                #     print(poly1.name, conf1, poly2.name, conf2)
                #     # obtain all properties
                #     # check N-N distance of poly1-conf > poly2-conf
                #     NN_dist1 = np.linalg.norm(PROP1['NN_v'])
                #     NN_dist2 = np.linalg.norm(PROP2['NN_v'])
                #     # determine angles made by NN_v and NN-BC_v
                #     # check that the pairs sum to 180
                #     p1_angle1 = PROP1['NN_BCN_1']
                #     p1_angle2 = PROP1['NN_BCN_2']
                #     p2_angle1 = PROP2['NN_BCN_1']
                #     p2_angle2 = PROP2['NN_BCN_2']
                #     # now check that the length of the long
                #     # vector and the short vector are commensurate
                #     # with an ideal trapezoid with the given angles
                #     # i.e. the extender vector determined by the
                #     # difference of the two NN_v (LHS) matches what is
                #     # expected by trig (RHS)
                #     extender_V_LHS = (NN_dist1 - NN_dist2) / 2
                #     test_angle = np.radians(180 - p2_angle1)
                #     extender_V_RHS = vector_length * np.cos(test_angle)
                #     tol = 2 * vector_std * np.cos(test_angle)
                #     # temporary for plotting
                #     NN_dist2s.append(NN_dist2)
                #     angles1_sums.append(p1_angle1 + p2_angle1)
                #     angles2_sums.append(p1_angle2 + p2_angle2)
                #     LHSs.append(extender_V_LHS)
                #     RHSs.append(extender_V_RHS)
                #     # run checks
                #     if NN_dist1 > NN_dist2:
                #         if np.isclose(p1_angle1 + p2_angle1, 180, rtol=0, atol=angle_tol):
                #             if np.isclose(p1_angle2 + p2_angle2, 180, rtol=0, atol=angle_tol):
                #                 if np.isclose(extender_V_LHS, extender_V_RHS,
                #                               rtol=0, atol=tol):
                #                     print(poly1.name, conf1, poly2.name, conf2)
                #                     print('NN_dists:', NN_dist1, NN_dist2)
                #                     print('angles:', p1_angle1,
                #                           p1_angle2, p2_angle1, p2_angle2)
                #                     print('sum angles:',
                #                           p1_angle1 + p2_angle1,
                #                           p1_angle2 + p2_angle2)
                #                     print('LHS:', extender_V_LHS)
                #                     print('cos(pi-d):', test_angle)
                #                     print('RHS:', extender_V_RHS)
                #                     print('bond_tol:', tol)
                #                     print('passed')
                #                     input()
                #                     passed_pairs.append((i, j, conf1, conf2,
                #                                          p1_angle1, p1_angle2,
                #                                          p2_angle1, p2_angle2,
                #                                          NN_dist1, NN_dist2))
                #         all_pairs.append((i, j, conf1, conf2,
                #                           p1_angle1, p1_angle2,
                #                           p2_angle1, p2_angle2,
                #                           NN_dist1, NN_dist2))
                #     # make some temp plots
                #     fig, ax = plt.subplots(figsize=(8, 5))
                #     X_range = (5, 20)
                #     width = 0.5
                #     density = False
                #     Y = NN_dist2s
                #     xtitle = 'NN_dist2'
                #     ax.axvline(x=NN_dist1, c='k')
                #     X_bins = np.arange(X_range[0], X_range[1], width)
                #     hist, bin_edges = np.histogram(a=Y, bins=X_bins,
                #                                    density=density)
                #     ax.bar(bin_edges[:-1],
                #            hist,
                #            align='edge',
                #            alpha=1.0, width=width,
                #            color='firebrick',
                #            edgecolor='k')
                #     ax.tick_params(axis='both', which='major', labelsize=16)
                #     ax.set_xlabel(xtitle, fontsize=16)
                #     if density is False:
                #         ax.set_ylabel('count', fontsize=16)
                #     elif density is True:
                #         ax.set_ylabel('frequency', fontsize=16)
                #     ax.set_xlim(X_range)
                #     ax.set_title(poly1.name + '_' + str(conf1))
                #     fig.tight_layout()
                #     fig.savefig('NN_dists_' + poly1.name + '_' + str(conf1) + '_'+ poly2.name + '.pdf', dpi=720,
                #                 bbox_inches='tight')
                #     plt.close()
                #
                #     fig, ax = plt.subplots(figsize=(8, 5))
                #     X_range = (90, 270)
                #     width = 5
                #     density = False
                #     xtitle = 'sum of angles [degrees]'
                #     ax.axvline(x=180 - angle_tol, c='k')
                #     ax.axvline(x=180 + angle_tol, c='k')
                #
                #     Y = angles1_sums
                #     X_bins = np.arange(X_range[0], X_range[1], width)
                #     hist, bin_edges = np.histogram(a=Y, bins=X_bins,
                #                                    density=density)
                #     ax.bar(bin_edges[:-1],
                #            hist,
                #            align='edge',
                #            alpha=0.6, width=width,
                #            color='r',
                #            edgecolor='k', label='angle1 sum')
                #     Y = angles2_sums
                #     X_bins = np.arange(X_range[0], X_range[1], width)
                #     hist, bin_edges = np.histogram(a=Y, bins=X_bins,
                #                                    density=density)
                #     ax.bar(bin_edges[:-1],
                #            hist,
                #            align='edge',
                #            alpha=0.6, width=width,
                #            color='b',
                #            edgecolor='k', label='angle2 sum')
                #
                #     ax.tick_params(axis='both', which='major', labelsize=16)
                #     ax.set_xlabel(xtitle, fontsize=16)
                #     if density is False:
                #         ax.set_ylabel('count', fontsize=16)
                #     elif density is True:
                #         ax.set_ylabel('frequency', fontsize=16)
                #     ax.set_xlim(X_range)
                #     ax.set_title(poly1.name + '_' + str(conf1))
                #     ax.legend()
                #     fig.tight_layout()
                #     fig.savefig('angles_dists_' + poly1.name + '_' + str(conf1) + '_'+ poly2.name + '.pdf', dpi=720,
                #                 bbox_inches='tight')
                #     plt.close()
                #
                #     fig, ax = plt.subplots(figsize=(5, 5))
                #     X_range = (0, 5)
                #     X = LHSs
                #     Y = RHSs
                #     ax.scatter(X, Y, alpha=0.6, color='r', edgecolor='k')
                #     ax.plot(np.linspace(0, X_range[1], 10),
                #             np.linspace(0, X_range[1], 10), alpha=0.2, c='k')
                #     ax.tick_params(axis='both', which='major', labelsize=16)
                #     ax.set_xlabel('LHS', fontsize=16)
                #     ax.set_ylabel('RHS', fontsize=16)
                #     ax.set_xlim(X_range)
                #     ax.set_ylim(X_range)
                #     ax.set_title(poly1.name + '_' + str(conf1))
                #     fig.tight_layout()
                #     fig.savefig('LHSRHS_parity_' + poly1.name + '_' + str(conf1) + '_'+ poly2.name + '.pdf', dpi=720,
                #                 bbox_inches='tight')
                #     plt.close()


if __name__ == "__main__":
    main()
