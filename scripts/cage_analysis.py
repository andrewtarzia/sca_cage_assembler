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
import pandas as pd
import numpy as np
import os
import json

from utilities import (
    convert_symm_names,
    convert_lig_names_from_cage,
    get_plottables,
    start_at_0,
)

from plotting import (
    colors_i_like,
    heatmap,
    annotate_heatmap,
)


def plot_heatmap_X_vs_Y(
    data,
    xname,
    xlabel,
    symmetries,
    Y_name,
    ylabel,
    filename,
    experimentals,
    ylim=None,
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


def plot_4Dheatmap_X_vs_Y(
    ydata,
    yname,
    ylabel,
    xname,
    xlabel,
    x2name,
    x2label,
    symmetries,
    filename,
    experimentals,
    ylim=None,
):

    def size_scaling(x, x_ranges):
        max_size = 200
        val = (x/x_ranges[1]) * max_size
        return val

    # Define range for size scaling.
    allx2s = set([ydata[cageset][x2name] for cageset in ydata])
    x2range = (0, np.ceil(max(allx2s)))
    if ylim is None:
        allys = [
            ydata[cageset][yname][cage]
            for cageset in ydata
            for cage in ydata[cageset][yname]
            if ydata[cageset][yname][cage] is not None
        ]
        yrange = (min(allys), max(allys))
    else:
        yrange = ylim

    fig, ax = plt.subplots()
    x_range = np.arange(len(ydata.keys()))
    y_range = np.arange(len(symmetries))
    X, Y = np.meshgrid(x_range, y_range)
    X = X.flatten()
    Y = Y.flatten()
    Xnames = [list(ydata.keys())[i] for i in X]
    Ynames = [symmetries[i] for i in Y]
    Cs = [0 for i in range(len(Xnames))]
    X2s = [0 for i in range(len(Xnames))]
    for cagesetname in ydata:
        _ids = [
            i for i in range(len(Xnames))
            if Xnames[i] == cagesetname
        ]
        for _id in _ids:
            # Other cage set property.
            X2s[_id] = size_scaling(
                ydata[cagesetname][x2name],
                x_ranges=x2range,
            )

        # Iterate over all cages in set.
        for cage in ydata[cagesetname][yname]:
            symm = cage.split('_')[-1]
            _id = [
                i for i in range(len(Xnames))
                if (Xnames[i] == cagesetname and Ynames[i] == symm)
            ][0]
            Cs[_id] = ydata[cagesetname][yname][cage]

            # If experimental, add a star!
            if cage in experimentals:
                ax.scatter(
                    X[_id]+0.25,
                    Y[_id]+0.25,
                    marker='*',
                    c='orange',
                    s=120,
                )

    p = ax.scatter(
        X,
        Y,
        c=Cs,
        s=X2s,
        # edgecolor='k',
        cmap='viridis',
        vmin=yrange[0],
        vmax=yrange[1],
    )

    # Add a colorbar
    cbar = fig.colorbar(p)
    cbar.ax.set_ylabel(ylabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(ydata.keys())))
    ax.set_yticks(np.arange(len(symmetries)))
    # ... and label them with the respective list entries.
    ax.set_xticklabels([
        f"{convert_lig_names_from_cage('_'.join(i.split('_')[1:]))}: "
        f"{round(ydata[i][xname], 2)}"
        for i in ydata
    ])
    ax.set_yticklabels([convert_symm_names(i) for i in symmetries])

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(
        top=True,
        bottom=False,
        labeltop=True,
        labelbottom=False
    )

    # Rotate the tick labels and set their alignment.
    plt.setp(
        ax.get_xticklabels(),
        rotation=-30,
        ha="right",
        rotation_mode="anchor"
    )
    ax.tick_params(which="minor", bottom=False, left=False)

    # Fake legend.
    for i in np.arange(x2range[0], x2range[1], step=0.5):
        if i == 0:
            continue
        ax.scatter(
            -100,
            -100,
            c='white',
            s=size_scaling(i, x2range),
            edgecolor='k',
            label=f'{i}'
        )
    ax.legend(fontsize=16)
    ax.set_xlim(min(X)-1, max(X)+1)
    ax.set_ylim(min(Y)-1, max(Y)+1)

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
         - maxMLlength: max Zn-N bond length in cage
         - formatione: xtb formation energy of cage
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
            'm_cube_shape': {
                i.name: cage_set.get_m_cube_shape(i.name)
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
            'maxMLlength': {
                i.name: cage_set.get_max_ML_distance(i.name)
                for i in cage_set.cages_to_build
            },
            'formatione': {
                i.name: cage_set.get_formation_energy(i.name)
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

    plottables = get_plottables(measures=measures, name=cage_set.name)

    _figure_path = 'figures'
    for p1, p2 in combinations(plottables, 2):
        p1_dict = plottables[p1]
        p1_filename = os.path.join(
            _figure_path, p1_dict['filename'],
        )
        p2_dict = plottables[p2]
        p2_filename = os.path.join(
            _figure_path, p2_dict['filename'],
        )
        p1_p2_filename = os.path.join(
            _figure_path, f'{cage_set.name}_{p1}_{p2}.pdf',
        )
        p2_p1_filename = os.path.join(
            _figure_path, f'{cage_set.name}_{p2}_{p1}.pdf'
        )
        pN1 = p1 if p2 == 'lsesum' else p2
        pN1_dict = p1_dict if pN1 == p1 else p2_dict
        pN2_dict = p1_dict if pN1 != p1 else p2_dict
        lse_filename = os.path.join(
            _figure_path,  f'{cage_set.name}_{pN1}_LSE.pdf'
        )

        if not exists(p1_filename):
            cage_set.plot_Y(
                data=p1_dict['data'],
                ylabel=p1_dict['ylabel'],
                ylim=p1_dict['ylim'],
                filename=p1_filename,
            )

        if not exists(p2_filename):
            cage_set.plot_Y(
                data=p2_dict['data'],
                ylabel=p2_dict['ylabel'],
                ylim=p2_dict['ylim'],
                filename=p2_filename,
            )

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


def write_csv(cage_sets, experimentals):
    """
    Write a .csv file with all numerical values for cage library.

    """

    cage_set_data = {}
    for cage_set in cage_sets:
        measures = cage_set_properties(cage_set)
        # For all HoCube sets, collect ligand data and
        # each cages measures.
        cage_set_data[cage_set.name] = {
            'AR': cage_set.ligand_aspect_ratio,
            'AB': cage_set.ligand_aspect_difference,
            'FAMM': cage_set.face_properties,
            'PPD': cage_set.ideal_pore_size,
            'LAR': cage_set.flex_properties['la_range'],
            'octop': measures['octop'],
            'rellsesum': start_at_0(data_dict=measures['lsesum']),
            'lsesum': measures['lsesum'],
            'minitors': measures['minitors'],
            'm_cube_shape': measures['m_cube_shape'],
            'maxcrplan': measures['maxcrplan'],
            'maxMLlength': measures['maxMLlength'],
            'porediam': measures['porediam'],
            'relformatione': start_at_0(
                data_dict=measures['formatione']
            ),
            'maxintangledev': measures['maxintangledev'],
        }

    tested_cage_sets = [
        key for key in cage_set_data
        if f'_{key}_' in '__'.join(experimentals)
    ]
    set_columns = [
        'cageset', 'symmetry', 'AR', 'AB',
        'FAMMi', 'FAMMii', 'FAMMiii', 'FAMMiv', 'FAMMv',
        'FAMMvi', 'FAMMvii',
        'PPD', 'LAR', 'octop',
        'lsesum', 'rellsesum', 'minitors', 'maxcrplan',
        'maxMLlength', 'porediam', 'relformatione',
        'maxintangledev', 'm_cube_shape', 'outcome', 'tested'
    ]
    symmetries = [
        'd2', 'th1', 'th2', 'td', 'tl', 's61', 's62', 'd31', 'd32',
        'd31n', 'd32n', 's41', 's42', 'c2v', 'c2h',
    ]
    dataframe = pd.DataFrame(columns=set_columns)

    for cagesetname in cage_set_data:
        csd = cage_set_data[cagesetname]
        for symm in symmetries:
            rowinfo = {i: None for i in set_columns}
            rowinfo['cageset'] = cagesetname
            rowinfo['symmetry'] = symm
            rowinfo['AR'] = csd['AR']
            rowinfo['AB'] = csd['AB']
            rowinfo['FAMMi'] = np.average(csd['FAMM']['i'])
            rowinfo['FAMMii'] = np.average(csd['FAMM']['ii'])
            rowinfo['FAMMiii'] = np.average(csd['FAMM']['iii'])
            rowinfo['FAMMiv'] = np.average(csd['FAMM']['iv'])
            rowinfo['FAMMv'] = np.average(csd['FAMM']['v'])
            rowinfo['FAMMvi'] = np.average(csd['FAMM']['vi'])
            rowinfo['FAMMvii'] = np.average(csd['FAMM']['vii'])
            rowinfo['PPD'] = csd['PPD']
            rowinfo['LAR'] = csd['LAR']
            cage_name = f'C_{cagesetname}_{symm}'
            if cage_name in experimentals:
                rowinfo['outcome'] = 1
            else:
                rowinfo['outcome'] = 0
            if cagesetname in tested_cage_sets:
                rowinfo['tested'] = 1
            else:
                rowinfo['tested'] = 0
            rowinfo['octop'] = csd['octop'][cage_name]
            rowinfo['lsesum'] = csd['lsesum'][cage_name]
            rowinfo['rellsesum'] = csd['rellsesum'][cage_name]
            rowinfo['minitors'] = csd['minitors'][cage_name]
            rowinfo['m_cube_shape'] = csd['m_cube_shape'][cage_name]
            rowinfo['maxcrplan'] = csd['maxcrplan'][cage_name]
            rowinfo['maxMLlength'] = csd['maxMLlength'][cage_name]
            rowinfo['porediam'] = csd['porediam'][cage_name]
            rowinfo['relformatione'] = csd['relformatione'][cage_name]
            rowinfo['maxintangledev'] = csd['maxintangledev'][cage_name]
            dataframe = dataframe.append(rowinfo, ignore_index=True)

    dataframe.to_csv('all_cage_csv_data.csv')


def write_xray_csv(xtal_cage_data):
    """
    Write a .csv file with all numerical values for cage library.

    """

    set_columns = [
        'cageset', 'symmetry', 'AR', 'AB', 'FAMM1', 'FAMM2', 'FAMM3', 'FAMM4',
        'FAMM5', 'PPD', 'LAR', 'octop',
        'lsesum', 'rellsesum', 'minitors', 'maxcrplan',
        'maxMLlength', 'porediam', 'relformatione',
        'maxintangledev', 'm_cube_shape', 'outcome', 'tested'
    ]
    dataframe = pd.DataFrame(columns=set_columns)

    for xtal in xtal_cage_data:
        csd = xtal_cage_data[xtal]

        rowinfo = {i: None for i in set_columns}
        rowinfo['cageset'] = xtal
        rowinfo['symmetry'] = None
        rowinfo['AR'] = None
        rowinfo['AB'] = None
        rowinfo['FAMM1'] = None
        rowinfo['FAMM2'] = None
        rowinfo['FAMM3'] = None
        rowinfo['FAMM4'] = None
        rowinfo['FAMM5'] = None
        rowinfo['PPD'] = None
        rowinfo['LAR'] = None
        rowinfo['outcome'] = 1
        rowinfo['tested'] = 1
        rowinfo['octop'] = csd['octop']
        rowinfo['lsesum'] = csd['lsesum']
        rowinfo['rellsesum'] = None
        rowinfo['minitors'] = csd['minitors']
        rowinfo['m_cube_shape'] = csd['m_cube_shape']
        rowinfo['maxcrplan'] = csd['maxcrplan']
        rowinfo['maxMLlength'] = csd['maxMLlength']
        rowinfo['porediam'] = csd['porediam']
        rowinfo['relformatione'] = None
        rowinfo['maxintangledev'] = csd['maxintangledev']
        dataframe = dataframe.append(rowinfo, ignore_index=True)

    dataframe.to_csv('all_xray_csv_data.csv')
