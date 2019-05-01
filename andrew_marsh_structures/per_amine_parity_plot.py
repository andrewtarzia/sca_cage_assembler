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


def main():
    """Run script.

    """
    if (not len(sys.argv) == 2):
        print("""
Usage: build_cages.py amine_type output_file wipe run_build
    output_file: file to output results
""")
        sys.exit()
    else:
        output_file = sys.argv[1]

    full_dataset = pd.read_csv(output_file)
    print(len(full_dataset))
    print(full_dataset.columns)
    list_of_aldes = sorted(list(set(full_dataset.bb1)))
    list_of_amines = sorted(list(set(full_dataset.bb2)))
    list_of_topos = sorted(list(set(full_dataset.topo)))
    print(list_of_amines, list_of_aldes, list_of_topos)
    print(len(list_of_amines))

    X_alde = 'aldehyde2'
    Y_alde = 'aldehyde3'
    property = {'axis_title': 'pore diameter [$\mathrm{\AA}$]',
                'column': 'p_diam_opt'}

    fig, ax = plt.subplots(figsize=(5, 5))
    for ami in list_of_amines:
        # print(ami)
        DF = full_dataset[full_dataset.bb2 == ami]
        # print(DF)
        # get topo DF
        for topo in list_of_topos:
            # print(topo)
            DFT = DF[DF.topo == topo]
            # print(DFT)
            # get X DF based on aldehyde
            X_DF = DFT[DFT.bb1 == X_alde]
            # print(X_DF)
            # get Y DF based on aldehyde
            Y_DF = DFT[DFT.bb1 == Y_alde]
            # print(Y_DF)
            X = float(X_DF[property['column']].iloc[0])
            Y = float(Y_DF[property['column']].iloc[0])
            print(X, Y)
            if bool(np.isclose(X, Y, rtol=0, atol=5)) is False:
                print(X_DF.bb1.iloc[0]+'_'+X_DF.bb2.iloc[0]+'_'+X_DF.topo.iloc[0])
                print(Y_DF.bb1.iloc[0]+'_'+Y_DF.bb2.iloc[0]+'_'+Y_DF.topo.iloc[0])
                input()
            ax.scatter(X, Y, c='mediumvioletred', alpha=0.8,
                       edgecolor='k', marker='o', s=80)

    lims = (0, round(max(full_dataset[property['column']]))+5)
    ax.plot(np.linspace(lims[0], lims[1], 3), np.linspace(lims[0], lims[1], 3),
            c='k', alpha=0.4)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(X_alde, fontsize=16)
    ax.set_ylabel(Y_alde, fontsize=16)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_title(property['axis_title'], fontsize=12)
    fig.tight_layout()
    fig.savefig("amine_parity_"+output_file.rstrip('.csv')+"_"+property['column']+".pdf",
                dpi=720, bbox_inches='tight')
    plt.close()
    sys.exit()


if __name__ == "__main__":
    main()
