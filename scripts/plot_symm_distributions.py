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

from utilities import convert_lig_names_from_cage, convert_symm_names


def yproperties():
    return {
        # 'octop': (
        #     r'min. $q_{\mathrm{oct}}$', (None, None), 'min'
        # ),
        'm_cube_shape': ('CU-8 cube measure', (-0.5, None), 'min'),
        'rellsesum': (
            r'rel. sum strain energy [kJ mol$^{-1}$]',
            (-10, 1000),
            'max'
        ),
        # 'lsesum': (
        #     r'sum strain energy [kJ mol$^{-1}$]',
        #     (None, None),
        #     'max'
        # ),
        # 'minitors': (
        #     r'min. imine torsion [degrees]', (None, None), 'min'
        # ),
        # 'maxcrplan': (
        #     r'max. ligand distortion [$\mathrm{\AA}$]',
        #     (None, None),
        #     'max'
        # ),
        # 'maxMLlength': (
        #     r'max. N-Zn bond length [$\mathrm{\AA}$]',
        #     (None, None),
        #     'max'
        # ),
        # 'porediam': (
        #     r'pore diameter [$\mathrm{\AA}$]', (None, None), 'min'
        # ),
        # 'relformatione': (
        #     r'rel. formation energy [kJ mol$^{-1}$]', (-10, 1000), 'max'
        # ),
        # 'maxintangledev': (
        #     r'max. interior angle deviation [degrees]',
        #     (None, None),
        #     'max'
        # ),
    }


def distribution_plot(df, col_name):
    _figure_path = 'figures'

    yprops = yproperties()

    symm_labels = convert_symm_names(no_symbol=True)

    fig, ax = plt.subplots(figsize=(8, 5))

    _x_positions = 0
    _x_names = []
    forms = []
    does_not_form = []
    for symm in symm_labels:

        if symm == 'tl':
            print('no')
            continue
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

    if col_name == 'm_cube_shape':
        ax.axhline(y=0, lw=2, linestyle='--', c='k')

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylabel(yprops[col_name][0], fontsize=16)
    ax.set_ylim(yprops[col_name][1])
    ax.set_xticks([i[0] for i in _x_names])
    ax.set_xticklabels([i[1] for i in _x_names], rotation=45)
    fig.legend(fontsize=16)
    fig.tight_layout()
    fig.savefig(
        os.path.join(_figure_path, f"sym_distribution_{col_name}.pdf"),
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def line_plot(df, col_name):
    _figure_path = 'figures'

    yprops = yproperties()
    symm_labels = convert_symm_names(no_symbol=True)

    cage_sets = {
        'cl1_quad2_5': '1',
        'cl1_quad2_16': '2',
        'cl1_quad2_12': '3',
        'cl1_quad2_3': '4',
        'cl1_quad2_8': '5',
        'cl1_quad2_2': '6',
    }

    feasible_syms = {
        'd2', 'th2', 'td', 'tl', 's62', 'd32', 'd31n', 'd32n'
    }

    for cs in cage_sets:
        fig, ax = plt.subplots(figsize=(8, 4))
        _x_positions = 0
        _x_names = []
        forms = []
        does_not_form = []
        unfeasible = []
        feasible = []

        cs_df = df[df['cageset'] == cs]

        for symm in symm_labels:
            if symm == 'tl':
                print('no')
                continue
            print_name = symm_labels[symm]
            set_df = cs_df[cs_df['symmetry'] == symm]
            _x_positions += 1
            _x_names.append((_x_positions, print_name))
            for i, row in set_df.iterrows():
                outcome = True if row['outcome'] == 1 else False
                y_val = row[col_name]

                if y_val is None:
                    continue
                if outcome:
                    forms.append((_x_positions, float(y_val)))
                else:
                    # does_not_form.append((_x_positions, float(y_val)))
                    if symm in feasible_syms:
                        feasible.append((_x_positions, float(y_val)))
                    else:
                        unfeasible.append((_x_positions, float(y_val)))

        # ax.scatter(
        #     x=[i[0] for i in does_not_form],
        #     y=[i[1] for i in does_not_form],
        #     c='gray',
        #     marker='o',
        #     alpha=1.0,
        #     s=120,
        #     label='does not form',
        # )
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
        ax.scatter(
            x=[i[0] for i in unfeasible],
            y=[i[1] for i in unfeasible],
            c='r',
            marker='D',
            edgecolors='k',
            alpha=1.0,
            s=120,
            label='unfeasible',
        )
        ax.scatter(
            x=[i[0] for i in feasible],
            y=[i[1] for i in feasible],
            c='gray',
            marker='P',
            edgecolors='k',
            alpha=1.0,
            s=120,
            label='feasible',
        )

        ax.axhline(y=0, lw=2, linestyle='--', c='k')

        # Set number of ticks for x-axis
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_ylabel(yprops[col_name][0], fontsize=16)
        ax.set_ylim(yprops[col_name][1])
        ax.set_xticks([i[0] for i in _x_names])
        ax.set_xticklabels([i[1] for i in _x_names], rotation=45)
        ax.set_title(f'{cage_sets[cs]}', fontsize=16)
        fig.legend(fontsize=16)
        fig.tight_layout()
        fig.savefig(
            os.path.join(_figure_path, f"sep_lsedist_{cs}.pdf"),
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
        # 'octop',
        'rellsesum',
        # 'minitors', 'lsesum',
        # 'maxcrplan', 'maxMLlength', 'porediam',
        # 'relformatione', 'maxintangledev',
        'm_cube_shape',
    ]
    for col_name in target_cols:
        distribution_plot(
            df=all_cage_properties,
            col_name=col_name,
        )
        if col_name == 'rellsesum':
            line_plot(
                df=all_cage_properties,
                col_name=col_name,
            )


if __name__ == "__main__":
    main()
