#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to calculate formation and reaction energies from GFN results.

Author: Andrew Tarzia

Date Created: 08 Mar 2019

"""

import os
import sys
import json
from collections import Counter
import matplotlib.pyplot as plt


def calculate_formation_energy(prod, react):
    '''Calculate formation energy of 'A' in a.u..

    Reaction formation energy == sum(product energy) - sum(reactant energy)

    Keyword arguments:
        prod (list) - list of product energies
        react (list) - list of reactant energies

    Returns:
        RFE (float) - Reaction formation energy in a.u.

    '''
    RFE = sum(prod) - sum(react)
    return RFE


def list_of_reactions():
    '''Helper function that organizes the list of reactions to calculate the
    energy of.

    Returns a list of dictionaries describing reactions.
    '''
    lor = [
    {'long-name': '(TFB)1(DACH)1(H2O)-1',
        'prod':['TFB1_DACH1_H2O-1', 'H2O', ],
        'react': ['TFB', 'DACH', ],
        'no.imine': 1,
        'size':'[1+1]'},
    {'long-name': '(TFB)1(DACH)2(H2O)-2',
        'prod':['TFB1_DACH2_H2O-2', 'H2O', 'H2O', ],
        'react': ['TFB', 'DACH', 'DACH', ],
        'no.imine': 2,
        'size':'[1+2]'},
    {'long-name': '(TFB)1(DACH)3(H2O)-3',
        'prod':['TFB1_DACH3_H2O-3', 'H2O', 'H2O', 'H2O', ],
        'react': ['TFB', 'DACH', 'DACH', 'DACH', ],
        'no.imine': 3,
        'size':'[1+3]'},
    {'long-name': '(TFB)2(DACH)3(H2O)-5',
        'prod':['TFB2_DACH3_H2O-5', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', ],
        'react': ['TFB', 'TFB', 'DACH', 'DACH', 'DACH', ],
        'no.imine': 5,
        'size':'[2+3]'},
    {'long-name': '(TFB)2(DACH)3(H2O)-6',
        'prod':['TFB2_DACH3_H2O-6', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', ],
        'react': ['TFB', 'TFB', 'DACH', 'DACH', 'DACH', ],
        'no.imine': 6,
        'size':'[2+3]'},
    {'long-name': '(TFB)2(DACH)4(H2O)-6',
        'prod':['TFB2_DACH4_H2O-6', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', ],
        'react': ['TFB', 'TFB', 'DACH', 'DACH', 'DACH', 'DACH', ],
        'no.imine': 6,
        'size':'[2+4]'},
    {'long-name': '(TFB)3(DACH)5(H2O)-9',
        'prod':['TFB3_DACH5_H2O-9', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', ],
        'react': ['TFB', 'TFB', 'TFB', 'DACH', 'DACH', 'DACH', 'DACH', 'DACH', ],
        'no.imine': 9,
        'size':'[3+5]'},
    {'long-name': '(TFB)4(DACH)6(H2O)-12',
        'prod':['TFB4_DACH6_H2O-12', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', ],
        'react': ['TFB', 'TFB', 'TFB', 'TFB', 'DACH', 'DACH', 'DACH', 'DACH', 'DACH', 'DACH', ],
        'no.imine': 12,
        'size':'[4+6]'},
    {'long-name': '(TFB)1(PDA)1(H2O)-1',
        'prod':['TFB1_PDA1_H2O-1', 'H2O', ],
        'react': ['TFB', 'PDA', ],
        'no.imine': 1,
        'size':'[1+1]'},
    {'long-name': '(TFB)1(PDA)2(H2O)-2',
        'prod':['TFB1_PDA2_H2O-2', 'H2O', 'H2O', ],
        'react': ['TFB', 'PDA', 'PDA', ],
        'no.imine': 2,
        'size':'[1+2]'},
    {'long-name': '(TFB)1(PDA)3(H2O)-3',
        'prod':['TFB1_PDA3_H2O-3', 'H2O', 'H2O', 'H2O', ],
        'react': ['TFB', 'PDA', 'PDA', 'PDA', ],
        'no.imine': 3,
        'size':'[1+3]'},
    {'long-name': '(TFB)2(PDA)3(H2O)-5',
        'prod':['TFB2_PDA3_H2O-5', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', ],
        'react': ['TFB', 'TFB', 'PDA', 'PDA', 'PDA', ],
        'no.imine': 5,
        'size':'[2+3]'},
    {'long-name': '(TFB)2(PDA)3(H2O)-6',
        'prod':['TFB2_PDA3_H2O-6', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', ],
        'react': ['TFB', 'TFB', 'PDA', 'PDA', 'PDA', ],
        'no.imine': 6,
        'size':'[2+3]'},
    {'long-name': '(TFB)2(PDA)4(H2O)-6',
        'prod':['TFB2_PDA4_H2O-6', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', ],
        'react': ['TFB', 'TFB', 'PDA', 'PDA', 'PDA', 'PDA', ],
        'no.imine': 6,
        'size':'[2+4]'},
    {'long-name': '(TFB)3(PDA)6(H2O)-9',
        'prod':['TFB3_PDA6_H2O-9', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', ],
        'react': ['TFB', 'TFB', 'TFB', 'PDA', 'PDA', 'PDA', 'PDA', 'PDA', 'PDA', ],
        'no.imine': 9,
        'size':'[3+6]'},
    {'long-name': '(TFB)4(PDA)6(H2O)-12',
        'prod':['TFB4_PDA6_H2O-12', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', ],
        'react': ['TFB', 'TFB', 'TFB', 'TFB', 'PDA', 'PDA', 'PDA', 'PDA', 'PDA', 'PDA', ],
        'no.imine': 12,
        'size':'[4+6]'},
    {'long-name': '(TFB)2(IPA)1(DACH)4(H2O)-8',
        'prod':['IC1', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', ],
        'react': ['TFB', 'TFB', 'DACH', 'DACH', 'DACH', 'DACH', 'IPA'],
        'no.imine': 8,
        'size':'[3+4]'},
    {'long-name': '(TFB)2(IPA)2(DACH)5(H2O)-10',
        'prod':['IC2', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', ],
        'react': ['TFB', 'TFB', 'DACH', 'DACH', 'DACH', 'DACH', 'DACH',  'IPA', 'IPA'],
        'no.imine': 10,
        'size':'[4+5]'},
    {'long-name': '(TFB)3(DACH)5(H2O)-9',
        'prod':['IC3', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', ],
        'react': ['TFB', 'TFB', 'TFB', 'DACH', 'DACH', 'DACH', 'DACH', 'DACH', ],
        'no.imine': 9,
        'size':'[3+5]'},
    {'long-name': '(TFB)3(IPA)1(DACH)5(H2O)-10',
        'prod':['IC4', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', ],
        'react': ['TFB', 'TFB', 'TFB', 'DACH', 'DACH', 'DACH', 'DACH', 'DACH', 'IPA'],
        'no.imine': 10,
        'size':'[4+5]'},
    {'long-name': '(TFB)3(IPA)1(DACH)6(H2O)-11',
        'prod':['IC5', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', ],
        'react': ['TFB', 'TFB', 'TFB', 'DACH', 'DACH', 'DACH', 'DACH', 'DACH', 'DACH', 'IPA'],
        'no.imine': 11,
        'size':'[4+6]'},
    {'long-name': '(TFB)3(IPA)1(DACH)7(H2O)-11',
        'prod':['IC6', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', ],
        'react': ['TFB', 'TFB', 'TFB', 'DACH', 'DACH', 'DACH', 'DACH', 'DACH', 'DACH', 'DACH', 'IPA'],
        'no.imine': 11,
        'size':'[4+7]'},
                    ]
    return lor


def flat_line(ax, x, y, w=0, C='k', m='x'):
    ax.plot([x-w, x, x+w], [y, y, y], c=C)
    ax.scatter(x, y, marker=m, c=C)


def figure_5(filename, RFEs):
    '''Recreate figure 5 in DOI: 10.1021/acs.chemmater.7b04323

    '''
    des_species = ['(TFB)1(DACH)1(H2O)-1', '(TFB)1(DACH)2(H2O)-2',
                   '(TFB)1(DACH)3(H2O)-3', '(TFB)2(DACH)3(H2O)-5',
                   '(TFB)2(DACH)3(H2O)-6', '(TFB)2(DACH)4(H2O)-6',
                   '(TFB)3(DACH)5(H2O)-9', '(TFB)4(DACH)6(H2O)-12']
    fig, ax = plt.subplots()
    # plot
    X_positions = {'[1+1]': 3, '[1+2]': 7, '[1+3]': 11, '[2+3]': 15,
                   '[2+4]': 19, '[3+5]': 23, '[4+6]': 27}
    X_values = []
    Y_values = []
    for i in RFEs:
        if i in des_species:
            X_values.append(X_positions[RFEs[i][2]])
            Y_values.append(RFEs[i][0] * 2625.50)
    for i, X in enumerate(X_values):
        flat_line(ax, x=X, y=Y_values[i], w=1.5, C='purple')
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('cluster size', fontsize=16)
    ax.set_ylabel('free energy of formation [kJ/mol]', fontsize=16)
    ax.set_xlim(0, 30)
    ax.set_ylim(-600, 20)
    ax.set_xticks(list(X_positions.values()))
    ax.set_xticklabels(list(X_positions.keys()))
    fig.tight_layout()
    fig.savefig(filename, dpi=720, bbox_inches='tight')


def figure_6(filename, RFEs):
    '''Recreate figure 6 in DOI: 10.1021/acs.chemmater.7b04323

    '''
    des_species = ['(TFB)1(PDA)1(H2O)-1', '(TFB)1(PDA)2(H2O)-2',
                   '(TFB)1(PDA)3(H2O)-3', '(TFB)2(PDA)3(H2O)-5',
                   '(TFB)2(PDA)3(H2O)-6', '(TFB)2(PDA)4(H2O)-6',
                   '(TFB)4(PDA)6(H2O)-12']
    fig, ax = plt.subplots()
    # plot
    X_positions = {'[1+1]': 3, '[1+2]': 7, '[1+3]': 11, '[2+3]': 15,
                   '[2+4]': 19, '[3+6]': 23, '[4+6]': 27}
    X_values = []
    Y_values = []
    for i in RFEs:
        if i in des_species:
            X_values.append(X_positions[RFEs[i][2]])
            Y_values.append(RFEs[i][0] * 2625.50)
    for i, X in enumerate(X_values):
        flat_line(ax, x=X, y=Y_values[i], w=1.5, C='orange')
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('cluster size', fontsize=16)
    ax.set_ylabel('free energy of formation [kJ/mol]', fontsize=16)
    ax.set_xlim(0, 30)
    ax.set_ylim(-600, 20)
    ax.set_xticks(list(X_positions.values()))
    ax.set_xticklabels(list(X_positions.keys()))
    fig.tight_layout()
    fig.savefig(filename, dpi=720, bbox_inches='tight')


def table_3(filename, RFEs):
    '''Recreate table 3 in DOI: 10.1021/acs.chemmater.7b04323 as a plot

    '''
    des_species = ['(TFB)2(IPA)1(DACH)4(H2O)-8',
                   '(TFB)2(IPA)2(DACH)5(H2O)-10',
                   '(TFB)3(DACH)5(H2O)-9',
                   '(TFB)3(IPA)1(DACH)5(H2O)-10',
                   '(TFB)3(IPA)1(DACH)6(H2O)-11',
                   '(TFB)3(IPA)1(DACH)7(H2O)-11',
                   '(TFB)4(DACH)6(H2O)-12']

    fig, ax = plt.subplots()
    # plot
    X_positions = {'IC1': 1, 'IC2': 2, 'IC3': 3, 'IC4': 4,
                   'IC5': 5, 'IC6': 6, 'CC3-R': 7}
    X_values = list(X_positions.values())
    Y_lit = [-45, -39.2, -40.1, -38.7, -39.8, -33.5, -44.1]
    ax.scatter(X_values, Y_lit, c='b', marker='o', edgecolor='k', s=80,
               label='literature')
    X_values = []
    Y_values = []
    for i in RFEs:
        if i in des_species:
            X_values.append(list(X_positions.values())[des_species.index(i)])
            Y_values.append(RFEs[i][0] * 2625.50 / RFEs[i][1])
    ax.scatter(X_values, Y_values, c='r', marker='X', edgecolor='k', s=80,
               label='GFN2-xTB')
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('structure', fontsize=16)
    ax.set_ylabel('free energy of formation \n per imine bondes [kJ/mol]',
                  fontsize=16)
    ax.set_xlim(0, 8)
    ax.set_ylim(-50, -30)
    ax.legend()
    ax.set_xticks(list(X_positions.values()))
    ax.set_xticklabels(list(X_positions.keys()))
    fig.tight_layout()
    fig.savefig(filename, dpi=720, bbox_inches='tight')


if __name__ == "__main__":
    if (not len(sys.argv) == 3):
        print("""
Usage: calculate_energies.py JSON energy
    JSON: JSON file to run analysis on.
    energy: which energy to use
        ('FE': free energy, 'TE': E_elec, 'HT': H(T=298K), 'GT': G(T=298K))
        ('RMSD': outputs the RMSD of all structures in LATEX table format)
""")
        sys.exit()
    else:
        JSON = sys.argv[1]
        EY = sys.argv[2]

    with open(JSON, 'r') as outfile:
        data = json.load(outfile)
    if EY != 'RMSD':
        new_JSON = JSON.replace('.json', '_RFE_'+EY+'.json')
        RFEs = {}
        # iterate over desired reactions and calculate RFE
        for R in list_of_reactions():
            # print(R)
            # print(R['prod'], R['react'])
            # print(Counter(R['prod']), Counter(R['react']))
            prod_energies = [float(data[i][EY]) for i in R['prod']]
            react_energies = [float(data[i][EY]) for i in R['react']]
            # print(prod_energies)
            # print(react_energies)
            RFE = calculate_formation_energy(prod_energies, react_energies)
            # print(RFE)
            RFEs[R['long-name']] = (RFE, R['no.imine'], R['size'])
            # input('ok?')
        for i in RFEs:
            print(i, RFEs[i][0])
        # for i in RFEs:
        #     print(i, RFEs[i][0] * 2625.50)
        with open(new_JSON, 'w') as outfile:
            json.dump(RFEs, outfile)
        fig5 = JSON.replace('.json', '_fig5_'+EY+'.pdf')
        fig6 = JSON.replace('.json', '_fig6_'+EY+'.pdf')
        tab3 = JSON.replace('.json', '_tab3_'+EY+'.pdf')
        # recreate plots from paper:
        # Figure 5
        figure_5(fig5, RFEs)
        # Figure 6
        figure_6(fig6, RFEs)
        # Table 3 as a plot
        table_3(tab3, RFEs)
    elif EY == 'RMSD':
        species = sorted(list(data.keys()))
        print(species)
        for i in species:
            print(i+' & '+str(round(data[i][EY], 4))+" \\\ ")