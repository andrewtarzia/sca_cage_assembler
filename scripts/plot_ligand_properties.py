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
from utilities import convert_lig_names_from_cage, read_lib


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


def adjust_lightness(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(
        c[0], max(0, min(1, amount * c[1])), c[2]
    )


def plot_MM_vs_AR(json_files, candms, expts, full=False, short=False, grouped=False):
    _figure_path = 'figures'

    if full:
        addon = '_full'
        ignored_faces = []
    else:
        ignored_faces = ['iii', 'v']
        addon = ''

    if short:
        dict_key = 'face_properties'
        suffix = '_short'
    else:
        dict_key = 'face_long_properties'
        suffix = ''

    if grouped:
        addon += '_grouped'
        ignored_faces = ['iii', 'iv', 'v', 'vi', 'vii']
        ncol = 2
        ylim = (0, 6)
    else:
        ncol = 5
        ylim = (0, 10)

    fig, ax = plt.subplots(figsize=(8, 5))
    stabs = {
        'i': {'ar': [], 'avgMDif': [], 'name': []},
        'ii': {'ar': [], 'avgMDif': [], 'name': []},
        'iii': {'ar': [], 'avgMDif': [], 'name': []},
        'iv': {'ar': [], 'avgMDif': [], 'name': []},
        'v': {'ar': [], 'avgMDif': [], 'name': []},
        'vi': {'ar': [], 'avgMDif': [], 'name': []},
        'vii': {'ar': [], 'avgMDif': [], 'name': []},
    }
    for i in json_files:
        cage_set = i.replace('_ligand_measures.json', '')
        if cage_set not in expts:
            continue
        print(f'doing {cage_set}')
        with open(i, 'r') as f:
            cage_set_data = json.load(f)
        lname = cage_set.replace('cl1_', '')
        for face in cage_set_data[dict_key]:
            stabs[face]['ar'].append(
                cage_set_data[dict_key][face][1]
            )
            stabs[face]['avgMDif'].append(
                cage_set_data[dict_key][face][0]
            )
            stabs[face]['name'].append(
                (
                    convert_lig_names_from_cage(lname, as_int=True),
                    convert_lig_names_from_cage(lname, as_sub=True)
                )
            )

    for face in stabs:
        if face in ignored_faces:
            continue
        c, m, lab = candms[face]

        XY1 = [
            (y, x, name) for y, x, name
            in sorted(
                zip(
                    stabs[face]['ar'],
                    stabs[face]['avgMDif'],
                    stabs[face]['name'],
                ),
                key=lambda pair: pair[0]
            )
        ]

        if face in ['i', 'ii']:
            ls = '-'
        else:
            ls = 'dashed'

        if grouped:
            groups = {
                'short120': {
                    'ids': [1, 2, 5],
                    'cnum': 1.2,
                    'labl': r'$\theta_{\mathrm{short}} = 120 \degree$',
                },
                'short60': {
                    'ids': [3, 4, 6],
                    'cnum': 0.8,
                    'labl': r'$\theta_{\mathrm{short}} = 60 \degree$',
                },
            }
            for group in groups:
                X = [i[0] for i in XY1 if i[2][0] in groups[group]['ids']]
                Y = [i[1] for i in XY1 if i[2][0] in groups[group]['ids']]
                name = [i[2] for i in XY1 if i[2][0] in groups[group]['ids']]
                labs = groups[group]['labl']
                ax.plot(
                    X,
                    Y,
                    c=adjust_lightness(c, groups[group]['cnum']),
                    # edgecolors='k',
                    marker='o',
                    alpha=1.0,
                    markersize=8,
                    linewidth=3,
                    linestyle='-',
                    label=f'{lab}: {labs}'
                )

                for x, y, s in zip(X, Y, name):
                    ax.text(s=s[1], x=x, y=y, fontsize=16)
            alpha = 0.4
            marker = None
        else:
            alpha = 1
            marker = m
            X = [i[0] for i in XY1]
            Y = [i[1] for i in XY1]
            name = [i[2][1] for i in XY1]
            ax.plot(
                X,
                Y,
                c=c,
                # edgecolors='k',
                marker=marker,
                alpha=alpha,
                markersize=8,
                linewidth=3,
                linestyle=ls,
                label=f'{lab}'
            )
            for x, y, s in zip(X, Y, name):
                ax.text(s=s, x=x, y=y, fontsize=16)

    # ax[0].tick_params(axis='both', which='major', labelsize=16)
    ax.tick_params(axis='both', which='major', labelsize=16)
    # ax[1].set_xlabel('aspect ratio [1:X]', fontsize=16)
    # ax.set_xlabel(
    #     r'aspect difference [$\mathrm{\AA}$]', fontsize=16
    # )
    ax.set_xlabel(
        (
            r'avg. $\Delta$N $\cdot\cdot\cdot$'
            r'N as AR [$\mathrm{\AA}$]'
        ),
        fontsize=16,
    )
    # ax[0].set_ylabel(r'avg. mismatch [$\mathrm{\AA}$]', fontsize=16)
    # ax.set_ylabel(r'mismatch [$\mathrm{\AA}$]', fontsize=16)
    ax.set_ylabel(
        (
            r'$\Delta$Zn$^{\mathrm{II}}\cdot\cdot\cdot$'
            r'Zn$^{\mathrm{II}}$ [$\mathrm{\AA}$]'
        ),
        fontsize=16,
    )
    # ax[0].set_xlim(0, 8)
    ax.set_xlim(0, 7)
    # ax[0].set_ylim(0, 7)
    ax.set_ylim(ylim)
    ax.legend(fontsize=16, ncol=ncol)
    fig.savefig(
        os.path.join(
            _figure_path, f'all_ligand_MM_vs_AR{addon}{suffix}.pdf'
        ),
        dpi=720,
        bbox_inches='tight'
    )

    plt.close()
    print(f'all_ligand_MM_vs_AR{addon}{suffix}.pdf')



def plot_MM_vs_NN(json_files, candms, expts, full=False, short=False, grouped=False):
    _figure_path = 'figures'

    if full:
        addon = '_full'
        ignored_faces = []
    else:
        ignored_faces = ['iii', 'v']
        addon = ''

    if short:
        dict_key = 'face_properties'
        suffix = '_short'
    else:
        dict_key = 'face_long_properties'
        suffix = ''

    if grouped:
        addon += '_grouped'
        ignored_faces = ['iii', 'iv', 'v', 'vi', 'vii']
        ncol = 2
        ylim = (0, 6)
    else:
        ncol = 5
        ylim = (0, 10)

    fig, ax = plt.subplots(figsize=(8, 5))
    stabs = {
        'i': {'nn': [], 'avgMDif': [], 'name': []},
        'ii': {'nn': [], 'avgMDif': [], 'name': []},
        'iii': {'nn': [], 'avgMDif': [], 'name': []},
        'iv': {'nn': [], 'avgMDif': [], 'name': []},
        'v': {'nn': [], 'avgMDif': [], 'name': []},
        'vi': {'nn': [], 'avgMDif': [], 'name': []},
        'vii': {'nn': [], 'avgMDif': [], 'name': []},
    }
    for i in json_files:
        cage_set = i.replace('_ligand_measures.json', '')
        if cage_set not in expts:
            continue
        print(f'doing {cage_set}')
        with open(i, 'r') as f:
            cage_set_data = json.load(f)
        lname = cage_set.replace('cl1_', '')
        for face in cage_set_data[dict_key]:
            stabs[face]['nn'].append(
                cage_set_data[dict_key][face][2]
            )
            stabs[face]['avgMDif'].append(
                cage_set_data[dict_key][face][0]
            )
            stabs[face]['name'].append(
                (
                    convert_lig_names_from_cage(lname, as_int=True),
                    convert_lig_names_from_cage(lname, as_sub=True)
                )
            )

    for face in stabs:
        if face in ignored_faces:
            continue
        c, m, lab = candms[face]

        XY1 = [
            (y, x, name) for y, x, name
            in sorted(
                zip(
                    stabs[face]['nn'],
                    stabs[face]['avgMDif'],
                    stabs[face]['name'],
                ),
                key=lambda pair: pair[0]
            )
        ]

        if face in ['i', 'ii']:
            ls = '-'
        else:
            ls = 'dashed'

        if grouped:
            groups = {
                'short120': {
                    'ids': [1, 2, 5],
                    'cnum': 1.2,
                    'labl': r'$\theta_{\mathrm{short}} = 120 \degree$',
                },
                'short60': {
                    'ids': [3, 4, 6],
                    'cnum': 0.8,
                    'labl': r'$\theta_{\mathrm{short}} = 60 \degree$',
                },
            }
            for group in groups:
                X = [i[0] for i in XY1 if i[2][0] in groups[group]['ids']]
                Y = [i[1] for i in XY1 if i[2][0] in groups[group]['ids']]
                name = [i[2] for i in XY1 if i[2][0] in groups[group]['ids']]
                labs = groups[group]['labl']
                ax.plot(
                    X,
                    Y,
                    c=adjust_lightness(c, groups[group]['cnum']),
                    # edgecolors='k',
                    marker='o',
                    alpha=1.0,
                    markersize=8,
                    linewidth=3,
                    linestyle='-',
                    label=f'{lab}: {labs}'
                )

                for x, y, s in zip(X, Y, name):
                    ax.text(s=s[1], x=x, y=y, fontsize=16)
            alpha = 0.4
            marker = None
        else:
            alpha = 1
            marker = m
            X = [i[0] for i in XY1]
            Y = [i[1] for i in XY1]
            name = [i[2][1] for i in XY1]
            ax.plot(
                X,
                Y,
                c=c,
                # edgecolors='k',
                marker=marker,
                alpha=alpha,
                markersize=8,
                linewidth=3,
                linestyle=ls,
                label=f'{lab}'
            )
            for x, y, s in zip(X, Y, name):
                ax.text(s=s, x=x, y=y, fontsize=16)

    # ax[0].tick_params(axis='both', which='major', labelsize=16)
    ax.tick_params(axis='both', which='major', labelsize=16)
    # ax[1].set_xlabel('aspect ratio [1:X]', fontsize=16)
    # ax.set_xlabel(
    #     r'aspect difference [$\mathrm{\AA}$]', fontsize=16
    # )
    ax.set_xlabel(
        (
            r'avg. $\Delta$N $\cdot\cdot\cdot$'
            r'N [$\mathrm{\AA}$]'
        ),
        fontsize=16,
    )
    # ax[0].set_ylabel(r'avg. mismatch [$\mathrm{\AA}$]', fontsize=16)
    # ax.set_ylabel(r'mismatch [$\mathrm{\AA}$]', fontsize=16)
    ax.set_ylabel(
        (
            r'$\Delta$Zn$^{\mathrm{II}}\cdot\cdot\cdot$'
            r'Zn$^{\mathrm{II}}$ [$\mathrm{\AA}$]'
        ),
        fontsize=16,
    )
    # ax[0].set_xlim(0, 8)
    ax.set_xlim(0, 7)
    # ax[0].set_ylim(0, 7)
    ax.set_ylim(ylim)
    ax.legend(fontsize=16, ncol=ncol)
    fig.savefig(
        os.path.join(
            _figure_path, f'all_ligand_MM_vs_NN{addon}{suffix}.pdf'
        ),
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
    plot_MM_vs_AR(json_files, candms, experimental_results, short=True)
    plot_MM_vs_AR(json_files, candms, experimental_results, full=True)
    plot_MM_vs_AR(json_files, candms, experimental_results, grouped=True)

    plot_MM_vs_NN(json_files, candms, experimental_results)
    plot_MM_vs_NN(json_files, candms, experimental_results, grouped=True)


if __name__ == "__main__":
    main()
