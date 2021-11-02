#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot ligand AR, flex and preferred face.

Author: Andrew Tarzia

Date Created: 08 Nov 2020

"""

from glob import glob
import json
import sys
import matplotlib.pyplot as plt

from plotting import colors_i_like
from utilities import read_lib

def plot_all_ligand_properties(json_files, candms, expts):
    print(expts)
    fig, ax = plt.subplots(figsize=(8, 5))
    for i in json_files:
        cage_set = i.replace('_ligand_measures.json', '')
        print(f'doing {cage_set}')
        with open(i, 'r') as f:
            cage_set_data = json.load(f)

        preferred_face = min(
            cage_set_data['face_properties'],
            key=cage_set_data['face_properties'].get
        )
        x = cage_set_data['ligand_aspect_ratio']
        y = cage_set_data['flex_properties']['la_range']
        c, m = candms[preferred_face]
        if cage_set in expts:
            exptl = expts[cage_set]
            if exptl == preferred_face:
                tc = 'k'
            else:
                tc = 'r'
            ax.text(x+0.03, y, exptl, c=tc, fontsize=16)
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
            label=f'{cm}'
        )

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('aspect ratio [1:X]', fontsize=16)
    ax.set_ylabel(
        r'long axis deviation [$\mathrm{\AA}$]',
        fontsize=16
    )
    ax.set_xlim(1.0, 2.1)
    ax.set_ylim(0.1, 2.1)
    ax.legend(fontsize=16, ncol=2)
    fig.savefig(
        'all_ligand_properties.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_MM_vs_AR(json_files, candms, expts):
    fig, ax = plt.subplots(figsize=(8, 5))
    stabs = {
        '1': {'ar': [], 'stab': []},
        '2': {'ar': [], 'stab': []},
        '3': {'ar': [], 'stab': []},
        '4': {'ar': [], 'stab': []},
        '5': {'ar': [], 'stab': []},
        '6': {'ar': [], 'stab': []},
        '7': {'ar': [], 'stab': []},
    }
    for i in json_files:
        cage_set = i.replace('_ligand_measures.json', '')
        print(f'doing {cage_set}')
        with open(i, 'r') as f:
            cage_set_data = json.load(f)
        for face in cage_set_data['face_properties']:
            stabs[face]['ar'].append(
                cage_set_data['ligand_aspect_ratio']
            )
            stabs[face]['stab'].append(
                cage_set_data['face_properties'][face]
            )

    for face in stabs:
        if face in ['1', '4', '6', '7']:
            continue
        c, m = candms[face]
        XY = [
            (y, x) for y, x
            in sorted(
                zip(stabs[face]['ar'], stabs[face]['stab']),
                key=lambda pair: pair[0]
            )
        ]
        print(XY)
        # input()
        X = [i[0] for i in XY]
        Y = [i[1] for i in XY]
        ax.plot(
            X,
            Y,
            c=c,
            # edgecolors='k',
            marker=m,
            alpha=1.0,
            markersize=10,
            linestyle='dashed',
            label=f'{face}'
        )

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('aspect ratio [1:X]', fontsize=16)
    ax.set_ylabel('avg. side mismatch [%]', fontsize=16)
    ax.set_xlim(1.0, 2.0)
    ax.set_ylim(0, None)
    ax.legend(fontsize=16, ncol=3)
    fig.savefig(
        'all_ligand_MM_vs_AR.pdf',
        dpi=720,
        bbox_inches='tight'
    )

    plt.close()


def main():
    first_line = (
        'Usage: plot_ligand_properties.py expt_lib_file'
    )
    if (not len(sys.argv) == 2):
        print(f"""
{first_line}

    expt_lib_file : (str)
        File containing experimental symmetry  information (XXXXX).

    """)
        sys.exit()
    else:
        expt_lib_file = sys.argv[1]

    json_files = glob('*_ligand_measures.json')
    expt_data = read_lib(expt_lib_file)

    # Dict of the cages with the face they form.
    experimental_results = {}
    for expt in expt_data:
        if expt_data[expt]['face'] is not None:
            experimental_results[expt] = expt_data[expt]['face']

    candms = {
        '1': (colors_i_like()[9], 'o'),
        '2': (colors_i_like()[4], 'X'),
        '3': (colors_i_like()[11], 's'),
        '4': (colors_i_like()[10], 'D'),
        '5': (colors_i_like()[3], 'P'),
        '6': (colors_i_like()[5], 'p'),
        '7': (colors_i_like()[7], '^'),
    }
    plot_all_ligand_properties(json_files, candms, experimental_results)
    plot_MM_vs_AR(json_files, candms, experimental_results)


if __name__ == "__main__":
    main()
