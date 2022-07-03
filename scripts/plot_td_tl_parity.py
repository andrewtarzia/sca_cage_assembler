#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot td vs tl parity over properties.

Author: Andrew Tarzia

Date Created: 23 Mar 2022

"""

import os
import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt

from utilities import convert_symm_names


def main():
    first_line = (
        'Usage: plot_td_tl_parity.py'
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

    int_angle = {}
    cu8 = {}
    lsesum = {}
    porediam = {}
    for i, row in all_cage_data.iterrows():
        symm = row['symmetry']
        if row['m_cube_shape'] is None:
            continue
        if row['maxintangledev'] is None:
            continue
        if symm not in ('tl', 'td'):
            continue
        cage_set = row['cageset']
        if cage_set not in int_angle:
            int_angle[cage_set] = {}
        if cage_set not in cu8:
            cu8[cage_set] = {}
        if cage_set not in lsesum:
            lsesum[cage_set] = {}
        if cage_set not in porediam:
            porediam[cage_set] = {}

        int_angle[cage_set][symm] = float(row['maxintangledev'])
        lsesum[cage_set][symm] = float(row['lsesum'])
        cu8[cage_set][symm] = float(row['m_cube_shape'])
        porediam[cage_set][symm] = float(row['porediam'])

    print('a', int_angle, '\n')
    print('b', lsesum, '\n')
    print('c', cu8, '\n')
    print('d', porediam, '\n')

    props = {
        'm_cube_shape': (
            cu8,
            'CU-8 cube measure',
            (0, 0.8),
        ),
        'lsesum': (
            lsesum,
            r'sum strain energy [kJ mol$^{-1}$]',
            (1500, 2200),
        ),
        'maxintangledev': (
            int_angle,
            r'max. interior angle deviation [degrees]',
            (0, 4),
        ),
        'porediam': (
            porediam,
            r'pore diameter [$\mathrm{\AA}$]',
            (0, 20),
        ),
    }

    for prop in props:
        fig, ax = plt.subplots(figsize=(5, 5))

        for cs in props[prop][0]:
            x = props[prop][0][cs]['tl']
            y = props[prop][0][cs]['td']
            ax.scatter(
                x,
                y,
                c='gold',
                edgecolors='k',
                marker='o',
                s=120,
                alpha=1.0,
            )

        fake_xs = np.linspace(props[prop][2][0], props[prop][2][1], 10)
        ax.plot(fake_xs, fake_xs, c='k', lw=2, linestyle='--')

        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel(convert_symm_names('tl'), fontsize=16)
        ax.set_ylabel(convert_symm_names('td'), fontsize=16)
        ax.set_title(props[prop][1], fontsize=16)
        ax.set_xlim(props[prop][2])
        ax.set_ylim(props[prop][2])

        fig.tight_layout()
        fig.savefig(
            os.path.join(_figure_path, f'td_tl_{prop}.pdf'),
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()


if __name__ == "__main__":
    main()
