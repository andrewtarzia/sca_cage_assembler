#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot ligand AR, flex and preferred face.

Author: Andrew Tarzia

Date Created: 08 Nov 2020

"""

import numpy as np
from glob import glob
import json
import sys
import matplotlib.pyplot as plt

from atools import colors_i_like
from stk.molecular.atoms.elements import Y


def main():
    first_line = (
        'Usage: plot_ligand_properties.py'
    )
    if (not len(sys.argv) == 1):
        print(f"""
{first_line}

    """)
        sys.exit()

    json_files = glob('*_ligand_measures.json')

    candms = {
        '1': (colors_i_like()[9], 'o'),
        '2': (colors_i_like()[10], 'X'),
        '3': (colors_i_like()[11], 's'),
        '4': (colors_i_like()[4], 'D'),
        '5': (colors_i_like()[3], 'P'),
    }
    fig, ax = plt.subplots(figsize=(8, 5))
    for i in json_files:
        cage_set = i.replace('_ligand_measures.json', '')
        print(f'doing {cage_set}')
        with open(i, 'r') as f:
            cage_set_data = json.load(f)

        print(cage_set_data)
        preferred_face = min(
            cage_set_data['face_properties'],
            key=cage_set_data['face_properties'].get
        )
        print(preferred_face)
        x = cage_set_data['ligand_aspect_ratio']
        y = cage_set_data['flex_properties']['la_range']
        c, m = candms[preferred_face]
        print(x, y, c, m)
        ax.scatter(
            x,
            y,
            c=c,
            edgecolors='k',
            marker=m,
            alpha=1.0,
            s=120
        )


    for cm in candms:
        c, m = candms[cm]
        ax.scatter(
            -100, -100,
            c=c,
            edgecolors='k',
            marker=m,
            alpha=1.0,
            s=120,
            label=f'face: {cm}'
        )

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('aspect ratio [1:X]', fontsize=16)
    ax.set_ylabel(
        r'long axis deviation [$\mathrm{\AA}$]',
        fontsize=16
    )
    ax.set_xlim(1.0, 3)
    ax.set_ylim(0.4, 2)
    ax.legend(fontsize=16)
    fig.savefig(
        'all_ligand_properties.pdf',
        dpi=720,
        bbox_inches='tight'
    )

    plt.close()


if __name__ == "__main__":
    main()
