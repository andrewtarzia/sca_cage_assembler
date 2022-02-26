#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot results of DFT calculations of ligand strain.

Author: Andrew Tarzia

Date Created: 15 Feb 2021


"""

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

from plotting import colors_i_like


def main():
    _figure_path = 'figures'

    csv_file = 'strain_energy_comparison.csv'
    data = pd.read_csv(csv_file)

    data['dft_kjpermol'] = data['dft'] * 2625.5

    fig, ax = plt.subplots(figsize=(5, 5))

    xray_data = data[data['from'] == 'crystal']
    calc_data = data[data['from'] == 'computation']
    print(xray_data)
    print(calc_data)

    ax.scatter(
        calc_data['xtb_kjpermol'],
        calc_data['dft_kjpermol'],
        c=colors_i_like()[10],
        edgecolors='k',
        marker='o',
        alpha=1.0,
        s=120,
        label='calculated',
    )

    ax.scatter(
        xray_data['xtb_kjpermol'],
        xray_data['dft_kjpermol'],
        c=colors_i_like()[11],
        edgecolors='k',
        marker='X',
        alpha=1.0,
        s=120,
        label='xray',
    )
    ax.plot(
        np.linspace(-10, 3000, 100), np.linspace(-10, 3000, 100),
        c='k'
    )
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(r'xTB strain energy [kJmol$^{-1}$]', fontsize=16)
    ax.set_ylabel(r'DFT strain energy [kJmol$^{-1}$]', fontsize=16)
    ax.set_xlim(0, 2500)
    ax.set_ylim(0, 2500)
    ax.legend(fontsize=16)
    fig.tight_layout()
    fig.savefig(
        os.path.join(_figure_path, 'strain_energy_comparison.pdf'),
        dpi=720, bbox_inches='tight'
    )

    ax.set_xlim(0, 600)
    ax.set_ylim(0, 600)
    ax.legend(fontsize=16)
    fig.tight_layout()
    fig.savefig(
        os.path.join(
            _figure_path, 'strain_energy_comparison_zoomed.pdf'
        ),
        dpi=720, bbox_inches='tight'
    )

    plt.close()


if __name__ == "__main__":
    main()
