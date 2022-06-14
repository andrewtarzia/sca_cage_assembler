#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot cube measure vs other properties.

Author: Andrew Tarzia

Date Created: 08 Nov 2020

"""

import os
import pandas as pd
import sys
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import numpy as np

from utilities import convert_symm_names


def main():
    first_line = (
        'Usage: plot_cube_vs_energy.py'
    )
    if (not len(sys.argv) == 1):
        print(f"""
{first_line}

""")
        sys.exit()
    else:
        pass

    _figure_path = 'figures'

    all_cage_data = pd.read_csv('all_cage_csv_data.csv')
    all_cage_data = all_cage_data.where(
        pd.notnull(all_cage_data), None
    )

    forms_x = []
    forms_ye = []
    forms_yangle = []
    no_forms_x = []
    no_forms_ye = []
    no_forms_yangle = []

    all_e_vs_x = {
        'd2': {'d': [], 'c': 'k', 'm': 'o'},
        # 'th1': {'d': [], 'c': 'r', 'm': 'D'},
        'th2': {'d': [], 'c': 'r', 'm': 'X'},
        # 'td': {'d': [], 'c': 'r', 'm': 'o'},
        'tl': {'d': [], 'c': 'gold', 'm': 'P'},
        # 's41': {'d': [], 'c': 'gold', 'm': 'X'},
        # 's42': {'d': [], 'c': 'gold', 'm': 'D'},
        # 's61': {'d': [], 'c': 'gray', 'm': 'X'},
        's62': {'d': [], 'c': 'gray', 'm': 'D'},
        # 'd31': {'d': [], 'c': 'skyblue', 'm': 'P'},
        'd32': {'d': [], 'c': 'skyblue', 'm': 'o'},
        # 'd31n': {'d': [], 'c': 'b', 'm': 'P'},
        # 'd32n': {'d': [], 'c': 'b', 'm': 'o'},
        # 'c2v': {'d': [], 'c': 'green', 'm': 'o'},
        # 'c2h': {'d': [], 'c': 'green', 'm': 'X'},
    }
    for i, row in all_cage_data.iterrows():
        if row['m_cube_shape'] is None:
            continue
        if row['maxintangledev'] is None:
            continue
        x = float(row['m_cube_shape'])
        ye = float(row['rellsesum'])
        yangle = float(row['maxintangledev'])

        if row['symmetry'] in all_e_vs_x:
            all_e_vs_x[row['symmetry']]['d'].append((x, ye))

        if int(row['outcome']) == 1:
            forms_x.append(x)
            forms_ye.append(ye)
            forms_yangle.append(yangle)
        elif int(row['outcome']) == 0:
            no_forms_x.append(x)
            no_forms_ye.append(ye)
            no_forms_yangle.append(yangle)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(
        no_forms_x,
        no_forms_ye,
        c='gray',
        edgecolors='none',
        marker='o',
        s=40,
        alpha=0.5,
        rasterized=True,
        label='does not form',
    )
    ax.scatter(
        forms_x,
        forms_ye,
        c='gold',
        edgecolors='k',
        marker='o',
        alpha=1.0,
        s=120,
        rasterized=True,
        label='forms',
    )

    # for x, y1, y2 in zip(no_forms_x, no_forms_yfe, no_forms_ye):

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('CU-8 cube measure', fontsize=16)
    ax.set_ylabel(
        r'rel. sum strain energy [kJmol$^{-1}$]',
        fontsize=16,
    )
    ax.set_xlim((-0.1, 2))
    # ax.set_ylim(yprops[col_name][1])
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(
        os.path.join(_figure_path, f"shape_vs_energies.pdf"),
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()

    fig, ax = plt.subplots(figsize=(8, 5))
    for symm in all_e_vs_x:
        points = np.array((
            [i[0] for i in all_e_vs_x[symm]['d']],
            [i[1] for i in all_e_vs_x[symm]['d']],
        )).T
        hull = ConvexHull(points)
        for i, simplex in enumerate(hull.simplices):
            if i == 0:
                label = convert_symm_names(symm)
            else:
                label = None
            ax.plot(
                points[simplex, 0],
                points[simplex, 1],
                c=all_e_vs_x[symm]['c'],
                lw=1,
                alpha=1.0,
                label=label,
            )

    ax.scatter(
        forms_x,
        forms_ye,
        c='none',
        edgecolors='k',
        marker='o',
        alpha=1.0,
        s=120,
        rasterized=True,
        label='forms',
    )

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('CU-8 cube measure', fontsize=16)
    ax.set_ylabel(
        r'rel. sum strain energy [kJmol$^{-1}$]',
        fontsize=16,
    )
    ax.set_xlim((-0.1, 2))
    # ax.set_ylim(yprops[col_name][1])
    ax.legend(fontsize=16, ncol=3)

    fig.tight_layout()
    fig.savefig(
        os.path.join(_figure_path, f"shape_vs_energies_col.pdf"),
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(
        no_forms_x,
        no_forms_yangle,
        c='gray',
        edgecolors='none',
        marker='o',
        s=40,
        alpha=0.5,
        rasterized=True,
        label='does not form',
    )
    ax.scatter(
        forms_x,
        forms_yangle,
        c='gold',
        edgecolors='k',
        marker='o',
        alpha=1.0,
        s=120,
        rasterized=True,
        label='forms',
    )

    # for x, y1, y2 in zip(no_forms_x, no_forms_yfe, no_forms_ye):

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('CU-8 cube measure', fontsize=16)
    ax.set_ylabel(
        r'max. interior angle deviation [degrees]', fontsize=16
    )
    ax.set_xlim(-0.1, 2)
    ax.set_ylim(-0.2, 6)
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(
        os.path.join(_figure_path, f"shape_vs_int_angle.pdf"),
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


if __name__ == "__main__":
    main()
