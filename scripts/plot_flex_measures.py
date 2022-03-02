#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot flex measures.

Author: Andrew Tarzia

Date Created: 12 Nov 2021

"""

from glob import glob
import json
import sys
import numpy as np
import matplotlib.pyplot as plt

from plotting import colors_i_like


def plot_dists(json_files, data_name, xlabel, xprops):
    num_mols = len(json_files)
    xwidth = xprops[2]
    xbins = np.arange(xprops[0], xprops[1], xwidth)
    fig, ax = plt.subplots(num_mols, 1, sharex=True, figsize=(8, 10))
    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=0)

    for i, _file in enumerate(json_files):
        lig_name = (
            _file.replace('_flex_measure.json', '').split('_')[1]
        )
        with open(_file, 'r') as f:
            data = json.load(f)
        measures = data[data_name]
        ax[i].text(xprops[0]+0.1, 0., lig_name, fontsize=16)
        ax[i].hist(
            x=measures,
            bins=xbins,
            density=True,
            # bottom=bottoms[i],
            histtype='stepfilled',
            # histtype='step',  fill=False,
            # stacked=True,
            linewidth=1.,
            facecolor='gold',
            alpha=1.0,
            # color=leg_info()[ami]['c'],
            # color='white',
            edgecolor='k',
            # label=str(lig_name),
        )
        ax[i].set_yticks([])
        # ax[i].set_ylim(-1, 1)

        ax[i].tick_params(left=False, bottom=False)
        # ax.set_ylim(0, ylim)

    ax[i].tick_params(labelsize=16, bottom=True)
    ax[i].set_xlabel(xlabel, fontsize=16)
    ax[i].set_ylabel('frequency', fontsize=16)
    fig.tight_layout()
    fig.savefig(
        f'flex_{data_name}_dists.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_ladev_planedev(json_files, series_props, series_map):
    fig, ax = plt.subplots(figsize=(8, 5))

    for i, _file in enumerate(json_files):
        with open(_file, 'r') as f:
            data = json.load(f)
        x = abs(
            max(data['long_axis_distances'])
            - min(data['long_axis_distances'])
        )
        y = abs(
            max(data['plane_deviations'])
            - min(data['plane_deviations'])
        )
        # measures = [i-min(measures) for i in measures]
        name = _file.replace('_flex_measure.json', '')
        prop = series_props[series_map[name]]
        ax.scatter(
            x=x,
            y=y,
            c=prop['c'],
            marker=prop['m'],
            alpha=1.0,
            edgecolor='k',
            s=160,
        )

    for prop in series_props:
        ax.scatter(
            [], [],
            c=series_props[prop]['c'],
            marker=series_props[prop]['m'],
            alpha=1.0,
            edgecolor='k',
            s=160,
            label=f'series {prop}',
        )

    ax.legend(fontsize=16)
    ax.tick_params(labelsize=16, bottom=True)
    ax.set_xlabel(r'long-axis deviation [$\mathrm{\AA}$]', fontsize=16)
    ax.set_ylabel(r'plane deviation [$\mathrm{\AA}$]', fontsize=16)
    fig.tight_layout()
    fig.savefig(
        'flex_comp.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_planedev_energy(json_files, series_props, series_map):
    num_mols = len(json_files)
    # xbins = np.arange(0, 15, xwidth)
    fig, ax = plt.subplots(
        num_mols, 1, sharex=True, sharey=True, figsize=(8, 10)
    )
    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=0)

    for i, _file in enumerate(json_files):
        lig_name = (
            _file.replace('_flex_measure.json', '').split('_')[1]
        )
        with open(_file, 'r') as f:
            data = json.load(f)
        x_measures = data['vector_angle']
        y_measures = [
            (i-min(data['energies']))*2625.5
            for i in data['energies']
        ]
        # measures = [i-min(measures) for i in measures]
        ax[i].text(0., 0., lig_name, fontsize=16)
        ax[i].scatter(
            x=x_measures,
            y=y_measures,
            facecolor='gold',
            alpha=1.0,
            edgecolor='k',
        )

        ax[i].tick_params(left=False, bottom=False)
    ax[i].tick_params(labelsize=16, bottom=True)
    ax[i].set_xlabel(
        'vecotr angle [degrees]', fontsize=16
    )
    ax[i].set_ylabel('rel. energy [a.u.]', fontsize=16)




    # # xbins = np.arange(0, 15, xwidth)
    # fig, ax = plt.subplots(figsize=(8, 5))
    # # Remove horizontal space between axes
    # fig.subplots_adjust(hspace=0)

    # for i, _file in enumerate(json_files):
    #     with open(_file, 'r') as f:
    #         data = json.load(f)
    #     x = abs(
    #         max(data['long_axis_distances'])
    #         - min(data['long_axis_distances'])
    #     )
    #     y = abs(
    #         max(data['plane_deviations'])
    #         - min(data['plane_deviations'])
    #     )
    #     # measures = [i-min(measures) for i in measures]
    #     ax.scatter(
    #         x=x,
    #         y=y,
    #         c='gold',
    #         alpha=1.0,
    #         edgecolor='k',
    #         s=160,
    #     )
    # ax.tick_params(labelsize=16, bottom=True)
    # ax.set_xlabel(r'long-axis deviation [$\mathrm{\AA}$]', fontsize=16)
    # ax.set_ylabel(r'plane deviation [$\mathrm{\AA}$]', fontsize=16)
    fig.tight_layout()
    fig.savefig(
        'flex_torsion_energy.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_ladev_angle(json_files, series_props, series_map):
    # xbins = np.arange(0, 15, xwidth)
    fig, ax = plt.subplots(figsize=(8, 5))

    for i, _file in enumerate(json_files):
        with open(_file, 'r') as f:
            data = json.load(f)
        x = abs(
            max(data['long_axis_distances'])
            - min(data['long_axis_distances'])
        )
        min_energy_angle = data['vector_angle'][0]
        y = abs(
            max(data['vector_angle'])
            - min(data['vector_angle'])
        )

        name = _file.replace('_flex_measure.json', '')
        prop = series_props[series_map[name]]
        ax.scatter(
            x=x,
            y=y,
            c=prop['c'],
            marker=prop['m'],
            alpha=1.0,
            edgecolor='k',
            s=160,
        )

    for prop in series_props:
        ax.scatter(
            [], [],
            c=series_props[prop]['c'],
            marker=series_props[prop]['m'],
            alpha=1.0,
            edgecolor='k',
            s=160,
            label=f'series {prop}',
        )

    ax.legend(fontsize=16)
    ax.tick_params(labelsize=16, bottom=True)
    ax.set_xlabel(r'long-axis deviation [$\mathrm{\AA}$]', fontsize=16)
    ax.set_ylabel(r'vector angle deviation [degrees]', fontsize=16)
    fig.tight_layout()
    fig.savefig(
        'flex_ladev_angle.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def main():
    first_line = (
        'Usage: plot_flex_measures.py'
    )
    if (not len(sys.argv) == 1):
        print(f"""{first_line}""")
        sys.exit()

    json_files = glob('*_flex_measure.json')
    series_map = {
        'quad2_1': '1',
        'quad2_12': '1',
        'quad2_2': '1',
        'quad2_3': '2',
        'quad2_8': '2',
        'quad2_9': '2',
        'quad2_10': '2',
        'quad2_5': '3',
        'quad2_16': '3',
        'quad2_17': '3',
    }
    series_props = {
        '1': {'c': colors_i_like()[9], 'm': 'o'},
        '2': {'c': colors_i_like()[4], 'm': 'X'},
        '3': {'c': colors_i_like()[11], 'm': 's'},
    }

    plot_dists(
        json_files=json_files,
        data_name='long_axis_distances',
        xlabel=r'long-axis deviation [$\mathrm{\AA}$]',
        xprops=(7.5, 25, 0.2),
    )
    plot_dists(
        json_files=json_files,
        data_name='vector_angle',
        xlabel=r'vector angle deviation [degrees]',
        xprops=(0., 90., 5.),
    )
    plot_ladev_planedev(json_files, series_props, series_map)
    plot_ladev_angle(json_files, series_props, series_map)
    plot_planedev_energy(json_files, series_props, series_map)


if __name__ == "__main__":
    main()
