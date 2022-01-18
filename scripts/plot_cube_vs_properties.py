#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot ligand strain energy vs formation energy for all cages.

Author: Andrew Tarzia

Date Created: 08 Nov 2020

"""

import numpy as np
import pandas as pd
import json
import sys
import matplotlib.pyplot as plt


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

    all_cage_data = pd.read_csv('all_cage_csv_data.csv')
    all_cage_data = all_cage_data.where(
        pd.notnull(all_cage_data), None
    )

    forms_x = []
    forms_yfe = []
    forms_ye = []
    forms_yangle = []
    no_forms_x = []
    no_forms_yfe = []
    no_forms_ye = []
    no_forms_yangle = []
    for i, row in all_cage_data.iterrows():
        if row['m_cube_shape'] is None:
            continue
        if row['maxintangledev'] is None:
            continue
        x = float(row['m_cube_shape'])
        yfe = float(row['m_cube_shape'])
        ye = float(row['rellsesum'])
        yangle = float(row['maxintangledev'])

        if int(row['outcome']) == 1:
            forms_x.append(x)
            forms_yfe.append(yfe)
            forms_ye.append(ye)
            forms_yangle.append(yangle)
        elif int(row['outcome']) == 0:
            no_forms_x.append(x)
            no_forms_yfe.append(yfe)
            no_forms_ye.append(ye)
            no_forms_yangle.append(yangle)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(
        no_forms_x,
        no_forms_ye,
        c='gray',
        edgecolors='none',
        marker='o',
        alpha=1.0,
        s=40,
        rasterized=True,
        label='does not form',
    )
    ax.scatter(
        forms_x,
        forms_ye,
        c='#4691C3',
        edgecolors='none',
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
    ax.set_ylabel('energies', fontsize=16)
    ax.set_xlim((-0.1, 2))
    # ax.set_ylim(yprops[col_name][1])
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(
        f"shape_vs_energies.pdf",
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
        alpha=1.0,
        s=40,
        rasterized=True,
        label='does not form',
    )
    ax.scatter(
        forms_x,
        forms_yangle,
        c='#4691C3',
        edgecolors='none',
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
    ax.set_ylim(-1, 10)
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(
        f"shape_vs_int_angle.pdf",
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


if __name__ == "__main__":
    main()
