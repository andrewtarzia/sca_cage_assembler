#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot parities of cage stability measures with xtals.

Author: Andrew Tarzia

"""

import pandas as pd
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
from utilities import read_lib, convert_lig_names_from_cage


def parity_plot(df, xray_df, col_name, pairings):
    _figure_path = 'figures'

    yprops = {
        'octop': (
            r'min. $q_{\mathrm{oct}}$', (0, 1), 'min'
        ),
        'm_cube_shape': ('CU-8 cube measure', (0, 1.6), 'min'),
        'rellsesum': (
            r'rel. sum strain energy [kJmol$^{-1}$]',
            (-10, 1000),
            'max'
        ),
        'lsesum': (
            r'sum strain energy [kJmol$^{-1}$]',
            (None, None),
            'max'
        ),
        'minitors': (
            r'min. imine torsion [degrees]', (160, 180), 'min'
        ),
        'maxcrplan': (
            r'max. ligand distortion [$\mathrm{\AA}$]',
            (50, 200),
            'max'
        ),
        'maxMLlength': (
            r'max. N-Zn bond length [$\mathrm{\AA}$]',
            (2.1, 2.4),
            'max'
        ),
        'porediam': (
            r'pore diameter [$\mathrm{\AA}$]', (5, 15), 'min'
        ),
        'relformatione': (
            r'rel. formation energy [kJmol$^{-1}$]', (-10, 1000), 'max'
        ),
        'maxintangledev': (
            r'max. interior angle deviation [degrees]',
            (0, 1),
            'max'
        ),
    }

    struct_map = {
        i: pairings[i]['xtal_struct_name']
        for i in pairings
    }

    fig, ax = plt.subplots(figsize=(5, 5))

    forms_df = df[df['outcome'] == 1]

    for i, row in forms_df.iterrows():
        cs = row['cageset']
        xray_name = struct_map[cs]
        xray_row = xray_df[xray_df['cageset'] == xray_name]
        calc_data = row[col_name]
        xray_data = xray_row[col_name]
        ax.scatter(
            calc_data,
            xray_data,
            c='#4691C3',
            edgecolors='k',
            marker='o',
            alpha=1.0,
            s=120,
        )
        ax.text(
            x=calc_data,
            y=xray_data,
            s=convert_lig_names_from_cage(cs[4:]),
            fontsize=16,
        )

    ax.plot(
        yprops[col_name][1], yprops[col_name][1],
        c='k',
        alpha=0.5,
        lw=2,
    )

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('calculated structure', fontsize=16)
    ax.set_ylabel('xray structure', fontsize=16)
    ax.set_title(yprops[col_name][0], fontsize=16)
    ax.set_xlim(yprops[col_name][1])
    ax.set_ylim(yprops[col_name][1])

    fig.tight_layout()
    fig.savefig(
        os.path.join(_figure_path, f"parities_{col_name}.pdf"),
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def main():
    first_line = (
        'Usage: plot_parities.py xray_data_path expt_lib_file'
    )
    if (not len(sys.argv) == 3):
        print(f"""
{first_line}

    xray_data_path : path to cray data csv
        ../xray_structures/analysis/all_xray_csv_data.csv

    expt_lib_file : (str)
        File containing experimental symmetry  information (XXXXX).

    """)
        sys.exit()
    else:
        xray_data_path = sys.argv[1]
        expt_lib_file = sys.argv[2]

    all_cage_properties = pd.read_csv('all_cage_csv_data.csv')
    all_xray_properties = pd.read_csv(xray_data_path)

    expt_data = read_lib(expt_lib_file)

    print(all_cage_properties.columns)
    target_cols = [
        'octop', 'minitors',
        'maxcrplan', 'maxMLlength', 'porediam',
        'maxintangledev',
        'm_cube_shape',
    ]
    for col_name in target_cols:
        parity_plot(
            df=all_cage_properties,
            xray_df=all_xray_properties,
            col_name=col_name,
            pairings=expt_data,
        )


if __name__ == "__main__":
    main()
