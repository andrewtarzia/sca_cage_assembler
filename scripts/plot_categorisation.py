#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot categorisation of cage stability measures.

Author: Andrew Tarzia

Date Created: 12 Feb 2021

"""

import pandas as pd
import os
import numpy as np
import sys
import matplotlib.pyplot as plt


def categorisation_plot(df, xray_df, col_name):
    _figure_path = 'figures'

    yprops = {
        'octop': (
            r'min. $q_{\mathrm{oct}}$', (None, None), 'min'
        ),
        'm_cube_shape': ('CU-8 cube measure', (-0.5, None), 'min'),
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
            r'min. imine torsion [degrees]', (None, None), 'min'
        ),
        'maxcrplan': (
            r'max. ligand distortion [$\mathrm{\AA}$]',
            (None, None),
            'max'
        ),
        'maxMLlength': (
            r'max. N-Zn bond length [$\mathrm{\AA}$]',
            (None, None),
            'max'
        ),
        'porediam': (
            r'pore diameter [$\mathrm{\AA}$]', (None, None), 'min'
        ),
        'relformatione': (
            r'rel. formation energy [kJmol$^{-1}$]', (-10, 1000), 'max'
        ),
        'maxintangledev': (
            r'max. interior angle deviation [degrees]',
            (None, None),
            'max'
        ),
    }

    dont_show = ['rellsesum', 'relformatione']

    dx = 0.15
    fig, ax = plt.subplots(figsize=(5, 5))

    forms_df = df[df['outcome'] == 1]
    does_not_form_df = df[df['outcome'] == 0]

    forms = list(forms_df[col_name])
    # names = list(forms_df['cageset'])
    does_not_form = list(does_not_form_df[col_name])
    xray = list(xray_df[col_name])
    for i in forms:
    # for i, j in zip(forms, names):
        ax.scatter(
            0.25+(dx*(np.random.random() - 0.5) * 2),
            i,
            c='#4691C3',
            edgecolors='k',
            marker='o',
            alpha=1.0,
            s=120,
        )
        # ax.text(
        #     x=0.25+(dx*(np.random.random() - 0.5) * 2),
        #     y=i,
        #     s=f'{j}',
        # )

    if col_name not in dont_show:
        for i in xray:
            ax.scatter(
                0.25+(dx*(np.random.random() - 0.5) * 2),
                i,
                c='#E74C3C',
                edgecolors='k',
                marker='X',
                alpha=1.0,
                s=120,
            )

    for i in does_not_form:
        ax.scatter(
            0.75+(dx*(np.random.random() - 0.5) * 2),
            i,
            c='#B8BEC3',
            edgecolors='k',
            marker='o',
            alpha=1.0,
            s=60,
        )

    if yprops[col_name][2] == 'max':
        if max(forms) < min(does_not_form):
            ax.axhline(y=max(forms), c='k', lw=2, linestyle='--')
    elif yprops[col_name][2] == 'min':
        if min(forms) > max(does_not_form):
            ax.axhline(y=min(forms), c='k', lw=2, linestyle='--')

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylabel(yprops[col_name][0], fontsize=16)
    ax.set_xlim((0, 1))
    ax.set_ylim(yprops[col_name][1])
    ax.set_xticks([0.25, 0.75])
    ax.set_xticklabels(['forms', 'does not form'])

    if col_name not in dont_show:
        ax.scatter(
            -100, 100,
            c='none',
            edgecolors='k',
            marker='o',
            alpha=1.0,
            s=120,
            label='calculated structure',
        )
        ax.scatter(
            -100, 100,
            c='#E74C3C',
            edgecolors='k',
            marker='X',
            alpha=1.0,
            s=120,
            label='xray structure',
        )
        ax.legend(fontsize=16)
    fig.tight_layout()
    fig.savefig(
        os.path.join(_figure_path, f"categorical_{col_name}.pdf"),
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def main():
    first_line = (
        'Usage: plot_categorisation.py xray_data_path'
    )
    if (not len(sys.argv) == 2):
        print(f"""
{first_line}

    xray_data_path : path to cray data csv
        ../xray_structures/analysis/all_xray_csv_data.csv

    """)
        sys.exit()
    else:
        xray_data_path = sys.argv[1]

    all_cage_properties = pd.read_csv('all_cage_csv_data.csv')
    all_xray_properties = pd.read_csv(xray_data_path)

    print(all_cage_properties.columns)
    target_cols = [
        'octop', 'rellsesum', 'minitors', 'lsesum',
        'maxcrplan', 'maxMLlength', 'porediam',
        'relformatione', 'maxintangledev',
        'm_cube_shape'
    ]
    for col_name in target_cols:
        categorisation_plot(
            all_cage_properties, all_xray_properties, col_name
        )


if __name__ == "__main__":
    main()
