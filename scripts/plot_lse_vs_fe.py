#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot ligand strain energy vs formation energy for all cages.

Author: Andrew Tarzia

Date Created: 08 Nov 2020

"""

from glob import glob
import json
import sys
import matplotlib.pyplot as plt


def main():
    first_line = (
        'Usage: plot_lse_vs_fe.py'
    )
    if (not len(sys.argv) == 1):
        print(f"""
{first_line}

    """)
        sys.exit()
    else:
        pass

    json_files = glob('*CS.json')

    lses = []
    form_eys = []
    for i in json_files:
        print(i)
        with open(i, 'r') as f:
            cage_data = json.load(f)
        print(cage_data.keys())
        for cage in cage_data:
            print(cage)
            fe = cage_data[cage]['fe_prop']
            li_data = cage_data[cage]['li_prop']
            lse_sum = sum([
                li_data['strain_energies'][i]
                for i in li_data['strain_energies']
            ])
            print(fe, li_data['strain_energies'], lse_sum)
            lses.append(lse_sum)
            form_eys.append(fe)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(
        lses,
        form_eys,
        c='k',
        edgecolors='k',
        marker='o',
        alpha=1.0,
        s=120
    )
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(r'sum strain energy [kJ/mol]', fontsize=16)
    ax.set_ylabel(r'formation energy [kJ/mol]', fontsize=16)
    # ax.set_xlim(1, i+3)
    # ax.set_ylim(ylim)

    fig.tight_layout()
    fig.savefig('lse_sum_vs_fe.pdf', dpi=720, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    main()
