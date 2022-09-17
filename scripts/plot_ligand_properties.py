#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot ligand AR, flex and preferred face.

Author: Andrew Tarzia

Date Created: 08 Nov 2020

"""

from glob import glob
import os
import json
import sys
import matplotlib.pyplot as plt

from plotting import colors_i_like
from utilities import read_lib


def plot_all_ligand_properties(json_files, candms, expts):
    _figure_path = 'figures'
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
        x = cage_set_data['ligand_aspect_difference']
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
    # ax.set_xlabel('aspect ratio [1:X]', fontsize=16)
    ax.set_xlabel(r'aspect difference [$\mathrm{\AA}$]', fontsize=16)
    ax.set_ylabel(
        r'long axis deviation [$\mathrm{\AA}$]',
        fontsize=16
    )
    ax.set_xlim(0., 8)
    ax.set_ylim(0.1, 2.1)
    ax.legend(fontsize=16, ncol=2)
    fig.savefig(
        os.path.join(_figure_path, 'all_ligand_properties.pdf'),
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_MM_vs_AR(json_files, candms, expts):
    _figure_path = 'figures'

    fig, ax = plt.subplots(figsize=(8, 5))
    stabs = {
        'i': {'ar': [], 'long_stab2': [], 'long_stab1': []},
        'ii': {'ar': [], 'long_stab2': [], 'long_stab1': []},
        'iii': {'ar': [], 'long_stab2': [], 'long_stab1': []},
        'iv': {'ar': [], 'long_stab2': [], 'long_stab1': []},
        'v': {'ar': [], 'long_stab2': [], 'long_stab1': []},
        'vi': {'ar': [], 'long_stab2': [], 'long_stab1': []},
        'vii': {'ar': [], 'long_stab2': [], 'long_stab1': []},
    }
    for i in json_files:
        cage_set = i.replace('_ligand_measures.json', '')
        if cage_set not in expts:
            continue
        print(f'doing {cage_set}')
        with open(i, 'r') as f:
            cage_set_data = json.load(f)
        for face in cage_set_data['face_properties']:
            stabs[face]['ar'].append(
                cage_set_data['ligand_aspect_difference']
            )
            stabs[face]['long_stab1'].append(
                cage_set_data['face_long_properties'][face][0]
            )
            stabs[face]['long_stab2'].append(
                cage_set_data['face_long_properties'][face][1]
            )

    for face in stabs:
        if face in ['iii', 'v']:
            continue
        c, m, lab = candms[face]
        XY1 = [
            (y, x) for y, x
            in sorted(
                zip(stabs[face]['ar'], stabs[face]['long_stab1']),
                key=lambda pair: pair[0]
            )
        ]
        XY2 = [
            (y, x) for y, x
            in sorted(
                zip(stabs[face]['ar'], stabs[face]['long_stab2']),
                key=lambda pair: pair[0]
            )
        ]
        X = [i[0] for i in XY1]
        Y = [i[1] for i in XY1]
        if face in ['i', 'ii']:
            ls = '-'
        else:
            ls = 'dashed'
        ax.plot(
            X,
            Y,
            c=c,
            # edgecolors='k',
            marker=m,
            alpha=1.0,
            markersize=8,
            linewidth=3,
            linestyle=ls,
            label=f'{lab}'
        )

    # ax[0].tick_params(axis='both', which='major', labelsize=16)
    ax.tick_params(axis='both', which='major', labelsize=16)
    # ax[1].set_xlabel('aspect ratio [1:X]', fontsize=16)
    ax.set_xlabel(
        r'aspect difference [$\mathrm{\AA}$]', fontsize=16
    )
    # ax[0].set_ylabel(r'avg. mismatch [$\mathrm{\AA}$]', fontsize=16)
    ax.set_ylabel(r'mismatch [$\mathrm{\AA}$]', fontsize=16)
    # ax[0].set_xlim(0, 8)
    ax.set_xlim(0, 7)
    # ax[0].set_ylim(0, 7)
    ax.set_ylim(0, 10)
    ax.legend(fontsize=16, ncol=5)
    fig.savefig(
        os.path.join(_figure_path, 'all_ligand_MM_vs_AR.pdf'),
        dpi=720,
        bbox_inches='tight'
    )

    plt.close()


def plot_MM_vs_AR_full(json_files, candms, expts):
    _figure_path = 'figures'

    fig, ax = plt.subplots(figsize=(8, 5))
    stabs = {
        'i': {'ar': [], 'long_stab2': [], 'long_stab1': []},
        'ii': {'ar': [], 'long_stab2': [], 'long_stab1': []},
        'iii': {'ar': [], 'long_stab2': [], 'long_stab1': []},
        'iv': {'ar': [], 'long_stab2': [], 'long_stab1': []},
        'v': {'ar': [], 'long_stab2': [], 'long_stab1': []},
        'vi': {'ar': [], 'long_stab2': [], 'long_stab1': []},
        'vii': {'ar': [], 'long_stab2': [], 'long_stab1': []},
    }
    for i in json_files:
        cage_set = i.replace('_ligand_measures.json', '')
        if cage_set not in expts:
            continue
        print(f'doing {cage_set}')
        with open(i, 'r') as f:
            cage_set_data = json.load(f)
        for face in cage_set_data['face_properties']:
            stabs[face]['ar'].append(
                cage_set_data['ligand_aspect_difference']
            )
            stabs[face]['long_stab1'].append(
                cage_set_data['face_long_properties'][face][0]
            )
            stabs[face]['long_stab2'].append(
                cage_set_data['face_long_properties'][face][1]
            )

    for face in stabs:
        c, m, lab = candms[face]
        XY1 = [
            (y, x) for y, x
            in sorted(
                zip(stabs[face]['ar'], stabs[face]['long_stab1']),
                key=lambda pair: pair[0]
            )
        ]
        XY2 = [
            (y, x) for y, x
            in sorted(
                zip(stabs[face]['ar'], stabs[face]['long_stab2']),
                key=lambda pair: pair[0]
            )
        ]

        X = [i[0] for i in XY1]
        Y = [i[1] for i in XY1]
        ls = '-'
        ax.plot(
            X,
            Y,
            c=c,
            # edgecolors='k',
            marker=m,
            alpha=1.0,
            markersize=8,
            linewidth=3,
            linestyle=ls,
            label=f'{face}-1'
        )

        X = [i[0] for i in XY2]
        Y = [i[1] for i in XY2]
        ls = 'dashed'
        ax.plot(
            X,
            Y,
            c=c,
            # edgecolors='k',
            marker=m,
            alpha=1.0,
            markersize=8,
            linewidth=3,
            linestyle=ls,
            label=f'{face}-2'
        )

    # ax[0].tick_params(axis='both', which='major', labelsize=16)
    ax.tick_params(axis='both', which='major', labelsize=16)
    # ax[1].set_xlabel('aspect ratio [1:X]', fontsize=16)
    ax.set_xlabel(
        r'aspect difference [$\mathrm{\AA}$]', fontsize=16
    )
    # ax[0].set_ylabel(r'avg. mismatch [$\mathrm{\AA}$]', fontsize=16)
    ax.set_ylabel(r'mismatch [$\mathrm{\AA}$]', fontsize=16)
    # ax[0].set_xlim(0, 8)
    ax.set_xlim(0, 7)
    # ax[0].set_ylim(0, 7)
    ax.set_ylim(0, 10)
    ax.legend(fontsize=16, ncol=5)
    fig.savefig(
        os.path.join(_figure_path, 'all_ligand_MM_vs_AR_full.pdf'),
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
            cs = expt_data[expt]['cage_set']
            experimental_results[cs] = expt_data[expt]['face']

    candms = {
        'i': (colors_i_like()[11], 'o', 'E1'),
        'ii': (colors_i_like()[4], 'X', 'E4'),
        'iii': (colors_i_like()[9], 's', 'E1'),
        'iv': (colors_i_like()[10], 'D', 'E5'),
        'v': (colors_i_like()[3], 'P', 'E1'),
        'vi': (colors_i_like()[5], 'p', 'E6'),
        'vii': (colors_i_like()[7], '^', 'E7'),
    }
    # plot_all_ligand_properties(
    #     json_files, candms, experimental_results,
    # )
    plot_MM_vs_AR(json_files, candms, experimental_results)
    plot_MM_vs_AR_full(json_files, candms, experimental_results)


if __name__ == "__main__":
    main()
