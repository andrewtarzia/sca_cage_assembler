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

from utilities import heatmap, annotate_heatmap, convert_symm_names

from atools import colors_i_like


def plot_heatmap_X_vs_Y(
    data,
    xname,
    xlabel,
    symmetries,
    Y_name,
    ylabel,
    ylim,
    filename,
    experimentals,
):

    col_names = {i: j[xname] for i, j in data.items()}
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

            # If experimental, add a star!
            if cage in experimentals:
                ax.scatter(
                    col-0.35,
                    row-0.4,
                    marker='*',
                    c='orange',
                    s=120,
                )

    im, cbar = heatmap(
        data=out_values,
        row_labels=[convert_symm_names(i) for i in symmetries],
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
    print(f'properties of: {cage_set.name}')
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


def analyse_cages(cage_sets, experimentals):

    AR_data = {}
    for cage_set in cage_sets:
        measures = cage_set_properties(cage_set)
        # For all HoCube sets, collect ligand data and
        # each cages measures.
        AR_data[cage_set.name] = {
            'AR': cage_set.ligand_aspect_ratio,
            'FAMM': cage_set.face_properties,
            'LAR': cage_set.flex_properties['la_range'],
            'octop': measures['octop'],
            'rellsesum': {
                i: (
                    measures['lsesum'][i]
                    - min(measures['lsesum'].values())
                )
                for i in measures['lsesum']
            },
            'minitors': measures['minitors'],
            'maxcrplan': measures['maxcrplan'],
            'maxdifffaceaniso': measures['maxdifffaceaniso'],
            'maxMLlength': measures['maxMLlength'],
            'relformatione': {
                i: (
                    measures['formatione'][i]
                    - min(measures['formatione'].values())
                )
                for i in measures['formatione']
            },
        }

    # Plot cage data as function of ligand data.
    tests = {
        # Test: (ylabel, ylim)
        'octop': (r'min. $q_{\mathrm{oct}}$', (0, 1)),
        'rellsesum': (r'rel. sum strain energy [kJ/mol]', (-4, 500)),
        'minitors': (r'min. imine torsion [degrees]', (0, 185)),
        'maxcrplan': (
            r'max. ligand distortion [$\mathrm{\AA}$]', (0, 185)
        ),
        'faceavgmismatches': (
            r'max. $\Delta$opposing face anisotropy [%]', (-10, 100)
        ),
        'maxMLlength': (
            r'max. N-Zn bond length [$\mathrm{\AA}$]', (2, 3)
        ),
        'relformatione': (
            r'rel. formation energy [kJ/mol]', (-10, 1000)
        ),
    }
    cs_tests = {
        # xlabel
        'AR': 'aspect ratio [1:X]',
        'LAR': r'long axis deviation [$\mathrm{\AA}$]',
        # 'FAMM': 'avg. side mismatch [%]',
    }
    if len(AR_data) > 0:
        for cs_test in cs_tests:
            # Plot vs cage properties.
            for t in tests:
                plot_heatmap_X_vs_Y(
                    data=AR_data,
                    xname=cs_test,
                    xlabel=cs_tests[cs_test],
                    symmetries=[
                        'o1', 'th1', 'th2', 't1', 's61',
                        's62', 'd31', 'd32', 'c2v', 'c2h',
                    ],
                    Y_name=t,
                    ylabel=tests[t][0],
                    ylim=tests[t][1],
                    filename=f'plotset_{t}_V{cs_test}.pdf',
                    experimentals=experimentals,
                )


def flat_line(ax, x, y, w=0, C='k', m='x', lw=2, label=None):
    ax.plot([x - w, x, x + w], [y, y, y], c=C, lw=lw)


def analyse_cage_sets(cage_sets):

    to_plot = {
        'AR': {
            'label': 'aspect ratio [1:X]',
        },
        'FAMM': {
            'label': 'avg. side mismatch [%]',
        },
        'LAR': {
            'label': r'long axis deviation [$\mathrm{\AA}$]',
        },
    }

    Cs = {
        '1': colors_i_like()[0],
        '2': colors_i_like()[1],
        '3': colors_i_like()[6],
        '4': colors_i_like()[3],
        '5': colors_i_like()[4],
    }
    M = 's'
    fig, ax = plt.subplots(figsize=(8, 5))
    for cs in cage_sets:
        ax.scatter(
            cs.ligand_aspect_ratio,
            cs.flex_properties['la_range'],
            c='none',
            edgecolors='k',
            marker=M,
            alpha=1.0,
            s=180
        )
        max_width = max(cs.face_properties.values())
        for fp in cs.face_properties:
            X = cs.ligand_aspect_ratio
            Y = cs.flex_properties['la_range']+(int(fp)-3)*0.02
            width = (cs.face_properties[fp]/max_width)*0.05
            flat_line(ax, x=X, y=Y, w=width, C=Cs[fp], lw=4)

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(to_plot['AR']['label'], fontsize=16)
    ax.set_ylabel(to_plot['LAR']['label'], fontsize=16)

    fig.tight_layout()
    fig.savefig(
        'cage_set_plot.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()
