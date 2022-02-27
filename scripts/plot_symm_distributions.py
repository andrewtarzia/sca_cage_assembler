#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot distributions of cage properties of each symm.

Author: Andrew Tarzia

Date Created: 12 Feb 2021

"""

import pandas as pd
import os
import sys
import matplotlib.pyplot as plt



def distribution_plot(df, col_name):
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

    symm_labels = {
        'd2': r'D$_2$',
        'th1': r'T$_{h, 1}$',
        'th2': r'T$_{h, 2}$',
        'td': r'T$_{\Delta}$',
        'tl': r'T$_{\Lambda}$',
        's41': r'S$_{4, 1}$',
        's42': r'S$_{4, 2}$',
        's61': r'S$_{6, 1}$',
        's62': r'S$_{6, 2}$',
        'd31': r'D$_{3, 1}$',
        'd32': r'D$_{3, 2}$',
        'd31n': r'D$_{3, 1n}$',
        'd32n': r'D$_{3, 2n}$',
        'c2v': r'C$_{2h}$',
        'c2h': r'C$_{2v}$',
    }

    fig, ax = plt.subplots(figsize=(8, 5))

    _x_positions = 0
    _x_names = []
    forms = []
    does_not_form = []
    for symm in symm_labels:
        print_name = symm_labels[symm]
        set_df = df[df['symmetry'] == symm]
        _x_positions += 1
        _x_names.append((_x_positions, print_name))
        cset_ys = []
        for i, row in set_df.iterrows():
            outcome = True if row['outcome'] == 1 else False
            y_val = row[col_name]
            if y_val is None:
                continue
            if outcome:
                forms.append((_x_positions, float(y_val)))
            else:
                does_not_form.append((_x_positions, float(y_val)))
            cset_ys.append(float(y_val))

        parts = ax.violinplot(
            cset_ys,
            [_x_positions],
            # points=200,
            vert=True,
            widths=0.8,
            showmeans=False,
            showextrema=False,
            showmedians=False,
            bw_method=0.5,
        )

        for pc in parts['bodies']:
            pc.set_facecolor('gray')
            pc.set_edgecolor('none')
            pc.set_alpha(0.3)


    ax.scatter(
        x=[i[0] for i in does_not_form],
        y=[i[1] for i in does_not_form],
        c='gray',
        marker='o',
        alpha=1.0,
        s=40,
        label='does not form',
    )
    ax.scatter(
        x=[i[0] for i in forms],
        y=[i[1] for i in forms],
        c='gold',
        edgecolors='k',
        marker='o',
        alpha=1.0,
        s=180,
        label='forms',
    )

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylabel(yprops[col_name][0], fontsize=16)
    ax.set_ylim(yprops[col_name][1])
    ax.set_xticks([i[0] for i in _x_names])
    ax.set_xticklabels([i[1] for i in _x_names])
    fig.legend(fontsize=16)
    fig.tight_layout()
    fig.savefig(
        os.path.join(_figure_path, f"sym_distribution_{col_name}.pdf"),
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def main():
    first_line = (
        'Usage: plot_symm_distributions.py'
    )
    if (not len(sys.argv) == 1):
        print(f"""
{first_line}
    """)
        sys.exit()
    else:
        pass

    all_cage_properties = pd.read_csv('all_cage_csv_data.csv')
    all_cage_properties = all_cage_properties.where(
        pd.notnull(all_cage_properties), None
    )

    target_cols = [
        'octop', 'rellsesum', 'minitors', 'lsesum',
        'maxcrplan', 'maxMLlength', 'porediam',
        'relformatione', 'maxintangledev',
        'm_cube_shape'
    ]
    for col_name in target_cols:
        distribution_plot(
            df=all_cage_properties,
            col_name=col_name,
        )


if __name__ == "__main__":
    main()
