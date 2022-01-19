#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot energies of cages in set.

Author: Andrew Tarzia

Date Created: 12 Feb 2021

"""

import pandas as pd
import os
import sys
import matplotlib.pyplot as plt


def plot_set_energies(data, filename):

    sets_to_plot = {
        'cl1_quad2_12': 'jd235',
        'cl1_quad2_8': 'jd257',
        'cl1_quad2_3': 'jd301',
        'cl1_quad2_16': 'jd326',
        'cl1_quad2_2': 'jd354',
        'cl1_quad2_5': 'jd370',
        'cl1_quad2_9': 'jd490',
    }

    fig, ax = plt.subplots(figsize=(8, 5))

    _x_positions = 0
    stab_energies = []
    _x_names = []
    for cageset in sets_to_plot:
        print_name = sets_to_plot[cageset]
        print(data[cageset])
        set_values = data[cageset]
        formed_symm = []
        all_energies = []
        for i in set_values:
            if set_values[i][1] is True:
                formed_symm.append((i, set_values[i]))
            all_energies.append(set_values[i][0])

        if len(formed_symm) != 1:
            raise ValueError(
                'Missing something; there should be one True case.'
            )
        print(formed_symm)
        print(all_energies)
        stabilisation_energy = formed_symm[0][1][0] - min(all_energies)
        print(stabilisation_energy)
        stab_energies.append(stabilisation_energy*2625.5)
        other_energies = [
            (i-min(all_energies))*2625.5 for i in all_energies
        ]
        _x_positions += 1
        _x_names.append((_x_positions, print_name))
        ax.scatter(
            x=[_x_positions for i in other_energies],
            y=other_energies,
            c='gray',
            s=40,
        )

    ax.scatter(
        x=range(1, _x_positions+1),
        y=stab_energies,
        c='r',
        s=180,
    )

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('ligand', fontsize=16)
    ax.set_ylabel(r'rel. energy [kJmol$^{-1}$]', fontsize=16)
    # ax.set_xlim((0, 1))
    ax.set_ylim(-0.1, None)
    ax.set_xticks([i[0] for i in _x_names])
    ax.set_xticklabels([i[1] for i in _x_names])

    fig.tight_layout()
    fig.savefig(
        filename,
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def main():
    first_line = (
        'Usage: plot_set_energies.py'
    )
    if (not len(sys.argv) == 1):
        print(f"""
{first_line}

    """)
        sys.exit()
    else:
        pass

    all_cage_properties = pd.read_csv('all_cage_csv_data.csv')

    # Define sets.
    sets = set(all_cage_properties['cageset'])

    # Get total energies from GFN.
    set_gfn_energies = {i: {} for i in sets}
    for i, row in all_cage_properties.iterrows():
        setname = row['cageset']
        symm = row['symmetry']
        forms = True if row['outcome'] == 1 else False
        ey_file = f'C_{setname}_{symm}_optc.ey'
        if not os.path.exists(ey_file):
            print(ey_file)
            continue
        with open(ey_file, 'r') as f:
            lines = f.readlines()
        total_energy = float(lines[0].rstrip())
        set_gfn_energies[setname][symm] = (total_energy, forms)

    # Get total energies from DFT.

    # Plot relative energy cf. more stable symmetry of expt symmetry.
    plot_set_energies(set_gfn_energies, 'set_energies_xtb.pdf')
    raise SystemExit('waiting on DFT')
    plot_set_energies(set_gfn_energies, 'set_energies_dft.pdf')


if __name__ == "__main__":
    main()
