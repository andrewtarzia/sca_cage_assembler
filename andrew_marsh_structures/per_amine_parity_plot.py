#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to compare the properties of cages from the same amine but two different
aldehydes.

Author: Andrew Tarzia

Date Created: 01 May 2019
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os.path import isfile
from stk import Cage, StructUnit2, StructUnit3
sys.path.insert(0, '/home/atarzia/thesource/')
from stk_functions import topo_2_property


def main():
    """Run script.

    """
    if (not len(sys.argv) == 3):
        print("""
Usage: per_amine_parity_plot.py output_file property
    output_file: file to output results
    property: which property you want to plot
        Defined:
            p_diam_opt - pore diameter in Angstrom
            bb_dist - distortion of building blocks
""")
        sys.exit()
    else:
        output_file = sys.argv[1]
        property_name = sys.argv[2]

    full_dataset = pd.read_csv(output_file)
    print(len(full_dataset))
    print(full_dataset.columns)
    list_of_aldes = sorted(list(set(full_dataset.bb1)))
    list_of_amines = sorted(list(set(full_dataset.bb2)))
    list_of_topos = sorted(list(set(full_dataset.topo)))
    list_of_names = sorted(list(full_dataset.name))
    # print(list_of_amines, list_of_aldes, list_of_topos)
    print(len(list_of_amines), 'x', len(list_of_aldes), 'x',
          len(list_of_topos), '-->', len(list_of_names))

    # define properties
    alde_pair_1 = ['aldehyde2', 'aldehyde3']
    alde_pair_2 = ['aldehyde1', 'aldehyde4']
    properties = {'p_diam_opt': {'axis_title': 'pore diameter [$\mathrm{\AA}$]',
                                 'column': 'p_diam_opt',
                                 'lims': (0,
                                          round(max(full_dataset['p_diam_opt']))+5)},
                  'bb_dist': {'axis_title': 'mean RMSD [$\mathrm{\AA}$]',
                              'column': None,  # needs to be calculated
                              'lims': (0, 6)}
                  }
    property = properties[property_name]
    if property_name == 'bb_dist':
        NEWCSV = output_file.replace('.csv', '_'+property_name+'.csv')
        # check if an output file exists.
        if isfile(NEWCSV):
            data_Frame = pd.read_csv(NEWCSV)
        else:
            # need to calculate the bb_distortion of all cages
            data_Frame = pd.DataFrame(columns=['NAME', 'value', 'topo', 'alde', 'amine'])
        done_list = list(set(data_Frame['NAME']))
        count = 1
        for name in list_of_names:
            if name in done_list:
                count += 1
                continue
            print('{} of {}'.format(count, len(list_of_names)))
            alde_name = name.split('_')[0]
            amine_name = '_'.join(name.split('_')[1:3])
            topo = name.split('_')[3]

            alde_dir = '/home/atarzia/projects/andrew_marsh_structures/precursor_lib/'
            alde_file = alde_dir+alde_name+'.mol'
            alde_struc = StructUnit3(alde_file, ['aldehyde'])
            amine_dir = '/data/atarzia/precursor_DBs/reaxys_sorted/'
            if '2f' in amine_name:
                amine_dir += 'amines2f/'
                amine_file = amine_dir+amine_name+'.mol'
                amine_struc = StructUnit2(amine_file, ['amine'])
            elif '3f' in amine_name:
                amine_dir += 'amines3f/'
                amine_file = amine_dir+amine_name+'.mol'
                amine_struc = StructUnit3(amine_file, ['amine'])

            topology = topo_2_property(topo, 'stk_func')
            cage = Cage([alde_struc, amine_struc], topology)
            cage.update_from_mol(name+'_opt.mol')
            data_Frame = data_Frame.append({'NAME': name,
                                            'value': cage.bb_distortion(),
                                            'amine': amine_name,
                                            'alde': alde_name,
                                            'topo': topo},
                                           ignore_index=True)
            count += 1
        # write dataFrame
        data_Frame.to_csv(NEWCSV,
                          index=False)

    fig, ax = plt.subplots(figsize=(5, 5))
    for ami in list_of_amines:
        if property['column'] is not None:
            DF = full_dataset[full_dataset.bb2 == ami]
            # get topo DF
            for topo in list_of_topos:
                DFT = DF[DF.topo == topo]
                # get X DF based on aldehyde
                X_DF = DFT[DFT.bb1 == alde_pair_1[0]]
                # get Y DF based on aldehyde
                Y_DF = DFT[DFT.bb1 == alde_pair_1[1]]
                X = float(X_DF[property['column']].iloc[0])
                Y = float(Y_DF[property['column']].iloc[0])
                ax.scatter(X, Y, c='mediumvioletred', alpha=0.6,
                           edgecolor='k', marker='o', s=80)
                # get X DF based on aldehyde
                X_DF = DFT[DFT.bb1 == alde_pair_2[0]]
                # get Y DF based on aldehyde
                Y_DF = DFT[DFT.bb1 == alde_pair_2[1]]
                X = float(X_DF[property['column']].iloc[0])
                Y = float(Y_DF[property['column']].iloc[0])
                ax.scatter(X, Y, c='orange', alpha=0.6,
                           edgecolor='k', marker='X', s=80)
        else:
            # plot from data_Frame
            DF = data_Frame[data_Frame['amine'] == ami]
            if DF.empty:
                continue
            for topo in list_of_topos:
                DFT = DF[DF.topo == topo]
                if DFT.empty:
                    continue
                # get X DF based on aldehyde
                X_DF = DFT[DFT.alde == alde_pair_1[0]]
                if X_DF.empty:
                    continue
                # get Y DF based on aldehyde
                Y_DF = DFT[DFT.alde == alde_pair_1[1]]
                if Y_DF.empty:
                    continue
                X = float(X_DF.value.iloc[0])
                Y = float(Y_DF.value.iloc[0])
                # print(X, Y)
                if abs(X-Y) > 1:
                    print(X_DF)
                    print(Y_DF)
                    input()
                ax.scatter(X, Y, c='mediumvioletred', alpha=0.6,
                           edgecolor='k', marker='o', s=80)
                # get X DF based on aldehyde
                X_DF = DFT[DFT.alde == alde_pair_2[0]]
                if X_DF.empty:
                    continue
                # get Y DF based on aldehyde
                Y_DF = DFT[DFT.alde == alde_pair_2[1]]
                if Y_DF.empty:
                    continue
                X = float(X_DF.value.iloc[0])
                Y = float(Y_DF.value.iloc[0])
                # print(X, Y)
                if abs(X-Y) > 1:
                    print(X_DF)
                    print(Y_DF)
                    input()
                ax.scatter(X, Y, c='orange', alpha=0.6,
                           edgecolor='k', marker='X', s=80)

    P = np.linspace(properties[property_name]['lims'][0],
                    properties[property_name]['lims'][1], 3)
    ax.plot(P, P, c='k', alpha=0.4)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('aldehyde 1/2', fontsize=16)
    ax.set_ylabel('aldehyde 3/4', fontsize=16)
    ax.set_xlim(properties[property_name]['lims'])
    ax.set_ylim(properties[property_name]['lims'])
    ax.set_title(property['axis_title'], fontsize=12)
    ax.scatter(X, Y, c='mediumvioletred', alpha=0.6,
               edgecolor='k', marker='o', s=80, label='aldehydes 2/3')
    ax.scatter(X, Y, c='orange', alpha=0.6,
               edgecolor='k', marker='X', s=80, label='aldehydes 1/4')
    ax.legend(fontsize=16)
    fig.tight_layout()
    fig.savefig("amine_parity_"+output_file.rstrip('.csv')+"_"+property_name+".pdf",
                dpi=720, bbox_inches='tight')
    plt.close()
    sys.exit()


if __name__ == "__main__":
    main()
