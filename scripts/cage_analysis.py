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
import json

from utilities import (
    heatmap,
    annotate_heatmap,
    convert_symm_names,
    convert_lig_names_from_cage,
    get_plottables,
    start_at_0,
)

from atools import colors_i_like


def plot_heatmap_X_vs_Y(
    data,
    xname,
    xlabel,
    symmetries,
    Y_name,
    ylabel,
    filename,
    experimentals,
):

    col_names = {
        i: (
            j[xname],
            convert_lig_names_from_cage('_'.join(i.split('_')[1:]))
        )
        for i, j in data.items()
    }
    col_names = {
        i: j for i, j in sorted(
            col_names.items(), key=lambda item: item[1][0]
        )
    }
    out_values = np.zeros((len(symmetries), len(col_names.values())))

    fig, ax = plt.subplots()

    na_points = []
    for name in data:
        col = list(col_names.keys()).index(name)
        # Iterate over all cages in set.
        for cage in data[name][Y_name]:
            Y = data[name][Y_name][cage]
            symm = cage.split('_')[-1]
            row = symmetries.index(symm)

            # Assign to zeros matrix.
            # Check if already assigned.
            if out_values[row][col] != 0:
                raise ValueError()

            if Y is None:
                na_points.append((row, col))
                # Set to zero and show marker.
                out_values[row][col] = np.nan
            else:
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
        col_labels=[
            f'{i[1]}: {round(i[0], 2)}' for i in col_names.values()
        ],
        ax=ax,
        cmap='Blues',
        cbarlabel=ylabel
    )
    annotate_heatmap(
        im=im,
        valfmt="{x:.1f}",
        fontsize=8,
        na_points=tuple(na_points)
    )

    ax.set_xlabel(xlabel)
    fig.tight_layout()
    fig.savefig(
        filename,
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()

    ax.set_xlabel(xlabel)
    fig.tight_layout()
    fig.savefig(
        filename,
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_X_vs_lse(
    xdata,
    ydata,
    xlabel,
    ylabel,
    filename,
    xlim=None,
    ylim=None,
):

    M = 'o'

    fig, ax = plt.subplots(figsize=(8, 5))

    xs = []
    ys = []
    for i, name in enumerate(xdata):
        xs.append(xdata[name])
        ys.append(ydata[name])
    ax.scatter(
        xs,
        ys,
        c=colors_i_like()[2],
        edgecolors='k',
        marker=M,
        alpha=1.0,
        s=180
    )

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

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
         - maxfacemetalpd:
         - maxintangledev: maximum deviation from 360 degrees
            (in degrees) of the interior angles of metal atoms in a
            face

    """

    cage_set.load_properties()

    # Get measures of all cages.
    if exists(cage_set.measures_file):
        with open(cage_set.measures_file, 'r') as f:
            measures = json.load(f)
    else:
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
            'porediam': {
                i.name: cage_set.get_pore_diameter(i.name)
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
            'maxfacemetalpd': {
                i.name: cage_set.get_max_face_metal_PD(i.name)
                for i in cage_set.cages_to_build
            },
            'maxintangledev': {
                i.name: cage_set.get_max_face_interior_angle_dev(
                    i.name
                )
                for i in cage_set.cages_to_build
            },
        }
        with open(cage_set.measures_file, 'w') as f:
            json.dump(measures, f, indent=4)

    print(f'properties of: {cage_set.name}')
    print('minimum OPs:', measures['octop'])
    print('sum LSE:', measures['lsesum'])
    print('min imine torsion:', measures['minitors'])
    print('max core planarities:', measures['maxcrplan'])
    print('max diff in face aniso:', measures['maxdifffaceaniso'])
    print('max metal-ligand distance:', measures['maxMLlength'])
    print('max face metal PD:', measures['maxfacemetalpd'])
    print('max face interior angle dev:', measures['maxintangledev'])
    print('pore diameters:', measures['porediam'])
    print('formation energies:', measures['formatione'])
    print('----------------------------------------------------')

    plottables = get_plottables(measures=measures, name=cage_set.name)

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

        if 'lsesum' in (p1, p2):
            pN1 = p1 if p2 == 'lsesum' else p2
            pN1_dict = p1_dict if pN1 == p1 else p2_dict
            pN2_dict = p1_dict if pN1 != p1 else p2_dict
            lse_filename = f'{cage_set.name}_{pN1}_LSE.pdf'
            if not exists(lse_filename):
                plot_X_vs_lse(
                    xdata=pN1_dict['data'],
                    ydata=pN2_dict['data'],
                    xlabel=pN1_dict['ylabel'],
                    ylabel=pN2_dict['ylabel'],
                    xlim=None,
                    ylim=None,
                    filename=lse_filename
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
            'PPD': cage_set.ideal_pore_size,
            'LAR': cage_set.flex_properties['la_range'],
            'octop': measures['octop'],
            'rellsesum': start_at_0(data_dict=measures['lsesum']),
            'minitors': measures['minitors'],
            'maxcrplan': measures['maxcrplan'],
            'maxdifffaceaniso': measures['maxdifffaceaniso'],
            'maxMLlength': measures['maxMLlength'],
            'porediam': measures['porediam'],
            'relformatione': start_at_0(
                data_dict=measures['formatione']
            ),
            'maxfacemetalpd': measures['maxfacemetalpd'],
            'maxintangledev': measures['maxintangledev'],
        }

    # Plot cage data as function of ligand data.
    tests = {
        # Test: ylabel
        'octop': r'min. $q_{\mathrm{oct}}$',
        'rellsesum': r'rel. sum strain energy [kJ/mol]',
        'minitors': r'min. imine torsion [degrees]',
        'maxcrplan': r'max. ligand distortion [$\mathrm{\AA}$]',
        'porediam': r'pore diameter [$\mathrm{\AA}$]',
        'maxdifffaceaniso': (
            r'max. $\Delta$opposing face anisotropy [%]'
        ),
        'maxfacemetalpd': (
            r'max. face metal planarity deviation [$\mathrm{\AA}$]'
        ),
        'maxintangledev': r'max. interior angle deviation [degrees]',
        'maxMLlength': r'max. N-Zn bond length [$\mathrm{\AA}$]',
        'relformatione': r'rel. formation energy [kJ/mol]',
    }
    cs_tests = {
        # xlabel
        'AR': 'aspect ratio [1:X]',
        'LAR': r'long axis deviation [$\mathrm{\AA}$]',
        'PPD': r'ideal pore size [$\mathrm{\AA}$]',
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
                    ylabel=tests[t],
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
