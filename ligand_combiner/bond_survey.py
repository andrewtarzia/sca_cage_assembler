#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to analyze the distribtuion of bonds in CSD subset.

Author: Andrew Tarzia

Date Created: 16 Apr 2019

"""

import sys
import pandas as pd
sys.path.insert(0, '/home/atarzia/thesource/')
import plotting


def main():
    if (not len(sys.argv) == 2):
        print("""
    Usage: bond_survey.py file
        file (str) - csv file to analyze
        """)
        sys.exit()
    else:
        file = sys.argv[1]

    data = pd.read_csv(file)

    fig, ax = plotting.histogram_plot_1(Y=data.DIST1, X_range=(1.8, 2.4),
                                        width=0.05, alpha=1.0, color='purple',
                                        edgecolor='k',
                                        xtitle='N-Pd bond distance [$\mathrm{\AA}$]')
    fig.tight_layout()
    fig.savefig('N_Pd_bond_survey.pdf', dpi=720,
                bbox_inches='tight')
    print('Mean = {} Angstrom'.format(data.DIST1.mean()))
    print('Standard Deviation = {} Angstrom'.format(data.DIST1.std()))


if __name__ == "__main__":
    main()
