#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot ligand strain energy vs formation energy for all cages.

Author: Andrew Tarzia

Date Created: 08 Nov 2020

"""

import numpy as np
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
        print(f'doing {i}')
        with open(i, 'r') as f:
            cage_set_data = json.load(f)

        print(
            f'there are {len(cage_set_data.keys())} cages in this set.'
        )
        for cage in cage_set_data:
            cage_data = cage_set_data[cage]
            print(
                f"doing cage: {cage}, optimized: "
                f"{cage_data['optimized']}."
            )
            if not cage_data['optimized']:
                continue
            fe = cage_data['fe_prop']
            li_data = cage_data['li_prop']
            lse_sum = sum([
                li_data['strain_energies'][i]
                for i in li_data['strain_energies']
            ])
            lses.append(lse_sum)
            form_eys.append(fe)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(
        lses,
        form_eys,
        c='gold',
        edgecolors='k',
        marker='o',
        alpha=1.0,
        s=80
    )
    ax.plot(
        lses,
        np.poly1d(np.polyfit(lses, form_eys, 1))(lses),
        c='k', lw=2
    )
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(r'sum strain energy [kJ/mol]', fontsize=16)
    ax.set_ylabel(r'formation energy [kJ/mol]', fontsize=16)

    fig.tight_layout()
    fig.savefig('lse_sum_vs_fe.pdf', dpi=720, bbox_inches='tight')

    ax.set_xlim(1550, 2500)
    ax.set_ylim(-2300, -1400)
    fig.savefig('lse_sum_vs_fe_z.pdf', dpi=720, bbox_inches='tight')

    plt.close()


if __name__ == "__main__":
    main()
