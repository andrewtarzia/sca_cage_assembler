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


def cage_set_properties(cage_set):
    """
    Collate cage-requiring properties of cage set.

    Analyses performed:
         - octop: min octahedral order parameter of zinc atoms
         - lsesum: sum of ligand xtb strain energies in cage
         - minitors: min imine torsion angle in cage
         - maxcrplan: max ligand core planarity in cage
         - maxdifffaceaniso: max percent face anisotropy of cage
         - maxMLlength: max Zn-N bond length in cage
         - formatione: xtb formation energy of cage

    """

    cage_set.load_properties()

    # Get measures of all cages.
    measures = {
        'octop': {
            i.name: cage_set.get_min_oct_op(i.name)
            for i in cage_set.cages_to_build
        },
        'lsesum': {
            i.name: cage_set.get_sum_lig_strain_energy(i.name)
            for i in cage_set.cages_to_build
        },
        'minitors': {
            i.name: cage_set.get_min_imine_torision(i.name)
            for i in cage_set.cages_to_build
        },
        'maxcrplan': {
            i.name: cage_set.get_max_core_planarity(i.name)
            for i in cage_set.cages_to_build
        },
        'maxdifffaceaniso': {
            i.name: cage_set.get_max_face_anisotropy(i.name)
            for i in cage_set.cages_to_build
        },
        'maxMLlength': {
            i.name: cage_set.get_max_ML_distance(i.name)
            for i in cage_set.cages_to_build
        },
        'formatione': {
            i.name: cage_set.get_formation_energy(i.name)
            for i in cage_set.cages_to_build
        },
    }
    print('minimum OPs:', measures['octop'])
    print('sum LSE:', measures['lsesum'])
    print('min imine torsion:', measures['minitors'])
    print('max core planarities:', measures['maxcrplan'])
    print('max diff in face aniso:', measures['maxdifffaceaniso'])
    print('max metal-ligand distance:', measures['maxMLlength'])
    print('formation energies:', measures['formatione'])
    print('----------------------------------------------------')

    plottables = {
        'octop': {
            'data': measures['octop'],
            'ylabel': r'min. $q_{\mathrm{oct}}$',
            'ylim': (0, 1),
            'filename': f'{cage_set.name}_minOPs.pdf'
        },
        'lsesum': {
            'data': {
                i: (
                    measures['lsesum'][i]
                    - min(measures['lsesum'].values())
                )
                for i in measures['lsesum']
            },
            'ylabel': r'rel. sum strain energy [kJ/mol]',
            'ylim': (-4, 500),
            'filename': f'{cage_set.name}_sumLSE.pdf'
        },
        'minitors': {
            'data': measures['minitors'],
            'ylabel': r'min. imine torsion [degrees]',
            'ylim': (0, 185),
            'filename': f'{cage_set.name}_mintors.pdf'
        },
        'maxcrplan': {
            'data': measures['maxcrplan'],
            'ylabel': r'max. core planarity [$\mathrm{\AA}$]',
            'ylim': (0, 185),
            'filename': f'{cage_set.name}_maxcrplane.pdf'
        },
        'maxdifffaceaniso': {
            'data': measures['maxdifffaceaniso'],
            'ylabel': r'max. $\Delta$opposing face anisotropy [%]',
            'ylim': (-10, 100),
            'filename': f'{cage_set.name}_maxfadiff.pdf'
        },
        'maxMLlength': {
            'data': measures['maxMLlength'],
            'ylabel': r'max. N-Zn bond length [$\mathrm{\AA}$]',
            'ylim': (2, 3),
            'filename': f'{cage_set.name}_maxmld.pdf'
        },
        'formatione': {
            'data': {
                i: (
                    measures['formatione'][i]
                    - min(measures['formatione'].values())
                )
                for i in measures['formatione']
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
