#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build HoCube library.

Author: Andrew Tarzia

Date Created: 27 Jan 2020

"""

from os.path import exists
from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np

from utilities import heatmap, annotate_heatmap


def plot_heatmap_X_vs_Y(
    data,
    xlabel,
    symmetries,
    Y_name,
    ylabel,
    ylim,
    filename
):

    col_names = {i: j['aspect_ratio'] for i, j in data.items()}
    col_names = {
        i: j for i, j in sorted(
            col_names.items(), key=lambda item: item[1]
        )
    }
    X_ = sorted(col_names.values())
    out_values = np.zeros((len(symmetries), len(X_)))

    fig, ax = plt.subplots()

    for name in data:
        col = list(col_names.keys()).index(name)
        # Iterate over all cages in set.
        for cage in data[name][Y_name]:
            Y = data[name][Y_name][cage]
            symm = cage.split('_')[-1]
            row = symmetries.index(symm)

            # Assign to zeros matrix.
            if out_values[row][col] != 0:
                raise ValueError()

            out_values[row][col] = Y

    im, cbar = heatmap(
        data=out_values,
        row_labels=symmetries,
        col_labels=[f'{round(i, 2)}' for i in col_names.values()],
        ax=ax,
        cmap='Blues',
        cbarlabel=ylabel
    )
    annotate_heatmap(im, valfmt="{x:.1f}", fontsize=8)
    ax.set_xlabel(xlabel)
    fig.tight_layout()
    fig.savefig(
        filename,
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def cage_set_analysis(cage_set):
    """
    Analyse cage set.

    Analysis performed:
         -
         -

    """

    cage_set.load_properties()
    built_prop = cage_set.built_cage_properties

    # Get measures of all cages.
    measures = {
        'oct_op': {},
        'lse_sum': {},
        'min_imine_torsions': {},
        'max_ligand_distortion': {},
        'max_diff_face_aniso': {},
        'max_ML_length': {},
        'formation_energies': {},
    }
    for C in cage_set.cages_to_build:
        C_data = built_prop[C.name]
        # Get minimium octahedral OP of the metal that is in the
        # complex building block.
        atom_no_of_interest = list(set([
            int(cage_set.complex_dicts[i]['metal_atom_no'])
            for i in cage_set.complex_dicts
        ]))
        C_OP = [
            C_data['op_prop'][str(i)]
            for i in atom_no_of_interest
        ]
        target_OPs = [
            i[j]['oct']
            for i in C_OP
            for j in i
        ]
        measures['oct_op'][C.name] = min(target_OPs)
        measures['lse_sum'][C.name] = sum([
            C_data['li_prop']['strain_energies'][i]
            for i in C_data['li_prop']['strain_energies']
        ])
        measures['min_imine_torsions'][C.name] = min([
            j
            for i in C_data['li_prop']['imine_torsions']
            for j in C_data['li_prop']['imine_torsions'][i]
        ])
        measures['max_ligand_distortion'][C.name] = max([
            C_data['li_prop']['core_planarities'][i]
            for i in C_data['li_prop']['core_planarities']
        ])
        measures['max_diff_face_aniso'][C.name] = max([
            100*((i[2] - i[3]) / i[2])
            for i in C_data['fa_prop']
        ])
        measures['max_ML_length'][C.name] = max([
            i for i in C_data['bl_prop']
        ])
        measures['formation_energies'][C.name] = C_data['fe_prop']

    print('minimum OPs:', measures['oct_op'])
    print('sum LSE:', measures['lse_sum'])
    print('min imine torsion:', measures['min_imine_torsions'])
    print('max core planarities:', measures['max_ligand_distortion'])
    print('max diff in face aniso:', measures['max_diff_face_aniso'])
    print('max metal-ligand distance:', measures['max_ML_length'])
    print('formation energies:', measures['formation_energies'])
    print('----------------------------------------------------')

    plottables = {
        'oct_op': {
            'data': measures['oct_op'],
            'ylabel': r'min. $q_{\mathrm{oct}}$',
            'ylim': (0, 1),
            'filename': f'{cage_set.name}_minOPs.pdf'
        },
        'lse_sum': {
            'data': {
                i: (
                    measures['lse_sum'][i]
                    - min(measures['lse_sum'].values())
                )
                for i in measures['lse_sum']
            },
            'ylabel': r'rel. sum strain energy [kJ/mol]',
            'ylim': (-4, 500),
            'filename': f'{cage_set.name}_sumLSE.pdf'
        },
        'min_imine_torsions': {
            'data': measures['min_imine_torsions'],
            'ylabel': r'min. imine torsion [degrees]',
            'ylim': (0, 185),
            'filename': f'{cage_set.name}_mintors.pdf'
        },
        'max_ligand_distortion': {
            'data': measures['max_ligand_distortion'],
            'ylabel': r'max. ligand distortion [$\mathrm{\AA}$]',
            'ylim': (0, 185),
            'filename': f'{cage_set.name}_maxdistortion.pdf'
        },
        'max_diff_face_aniso': {
            'data': measures['max_diff_face_aniso'],
            'ylabel': r'max. $\Delta$opposing face anisotropy [%]',
            'ylim': (-10, 100),
            'filename': f'{cage_set.name}_maxfadiff.pdf'
        },
        'max_ML_length': {
            'data': measures['max_ML_length'],
            'ylabel': r'max. N-Zn bond length [$\mathrm{\AA}$]',
            'ylim': (2, 3),
            'filename': f'{cage_set.name}_maxmld.pdf'
        },
        'relfe': {
            'data': {
                i: (
                    measures['formation_energies'][i]
                    - min(measures['formation_energies'].values())
                )
                for i in measures['formation_energies']
            },
            'ylabel': r'rel. formation energy [kJ/mol]',
            'ylim': (-10, 1000),
            'filename': f'{cage_set.name}_relfe.pdf'
        },
    }

    for p1, p2 in combinations(plottables, 2):
        p1_dict = plottables[p1]
        if not exists(p1_dict['filename']):
            cage_set.plot_Y(
                data=p1_dict['data'],
                ylabel=p1_dict['ylabel'],
                ylim=p1_dict['ylim'],
                filename=p1_dict['filename']
            )
        p2_dict = plottables[p2]
        if not exists(p2_dict['filename']):
            cage_set.plot_Y(
                data=p2_dict['data'],
                ylabel=p2_dict['ylabel'],
                ylim=p2_dict['ylim'],
                filename=p2_dict['filename']
            )

        p1_p2_filename = f'{cage_set.name}_{p1}_{p2}.pdf'
        if not exists(p1_p2_filename):
            cage_set.plot_Y_C(
                data=p1_dict['data'],
                ylabel=p1_dict['ylabel'],
                ylim=p1_dict['ylim'],
                data_C=p2_dict['data'],
                clabel=p2_dict['ylabel'],
                clim=p2_dict['ylim'],
                filename=p1_p2_filename
            )
        p2_p1_filename = f'{cage_set.name}_{p2}_{p1}.pdf'
        if not exists(p2_p1_filename):
            cage_set.plot_Y_C(
                data=p2_dict['data'],
                ylabel=p2_dict['ylabel'],
                ylim=p2_dict['ylim'],
                data_C=p1_dict['data'],
                clabel=p1_dict['ylabel'],
                clim=p1_dict['ylim'],
                filename=p2_p1_filename
            )

    return measures


def analyse_cages(cage_sets):

    AR_data = {}
    for cage_set in cage_sets:
        measures = cage_set_analysis(cage_set)
        # For all HoCube sets, collect ligand aspect ratio data.
        AR_data[cage_set.name] = {
            'min_OPs': measures['oct_op'],
            'lse_sum': {
                i: (
                    measures['lse_sum'][i] -
                    min(measures['lse_sum'].values())
                )
                for i in measures['lse_sum']
            },
            'min_tor': measures['min_imine_torsions'],
            'max_dis': measures['max_ligand_distortion'],
            'aspect_ratio': cage_set.ligand_aspect_ratio,
            'max_rfa': measures['max_diff_face_aniso'],
            'max_mld': measures['max_ML_length'],
            'rel_fe': {
                i: (
                    measures['formation_energies'][i]
                    - min(measures['formation_energies'].values())
                )
                for i in measures['formation_energies']
            },
        }

    # Plot ligand aspect ratio data.
    tests = {
        # Test: (ylabel, ylim)
        'min_OPs': (r'min. $q_{\mathrm{oct}}$', (0, 1)),
        'lse_sum': (r'rel. sum strain energy [kJ/mol]', (-4, 500)),
        'min_tor': (r'min. imine torsion [degrees]', (0, 185)),
        'max_dis': (
            r'max. ligand distortion [$\mathrm{\AA}$]', (0, 185)
        ),
        'max_rfa': (
            r'max. $\Delta$opposing face anisotropy [%]', (-10, 100)
        ),
        'max_mld': (
            r'max. N-Zn bond length [$\mathrm{\AA}$]', (2, 3)
        ),
        'rel_fe': (r'rel. formation energy [kJ/mol]', (-10, 1000)),
    }
    if len(AR_data) > 0:
        for t in tests:
            plot_heatmap_X_vs_Y(
                data=AR_data,
                xlabel='ligand aspect ratio [1:X]',
                symmetries=[
                    'o1', 'th1', 'th2', 't1', 's61', 's62', 'd31',
                    'd32', 'c2v', 'c2h',
                ],
                Y_name=t,
                ylabel=tests[t][0],
                ylim=tests[t][1],
                filename=f'plotset_{t}_VA.pdf'
            )
