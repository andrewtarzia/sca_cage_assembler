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


def plot_ladev_planedev(json_files):
    fig, ax = plt.subplots(figsize=(8, 5))

    prop = {'c': 'gold', 'm': 'o'}

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
        ax.scatter(
            x=x,
            y=y,
            c=prop['c'],
            marker=prop['m'],
            alpha=1.0,
            edgecolor='k',
            s=160,
        )

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


def plot_vectangle_energy_permol(json_files):

    _to_plot = ('2', '3', '5', '8', '12', '16')

    for i, _file in enumerate(json_files):
        lig_name = (
            _file.replace('_flex_measure.json', '').split('_')[1]
        )
        if lig_name not in _to_plot:
            continue

        fig, ax = plt.subplots(figsize=(8, 5))
        with open(_file, 'r') as f:
            data = json.load(f)
        x_measures = data['vector_angle']
        y_measures = [
            (i-min(data['energies']))*2625.5
            for i in data['energies']
        ]
        ax.scatter(
            x=x_measures,
            y=y_measures,
            facecolor='gold',
            alpha=1.0,
            edgecolor='k',
            s=120,
        )

        ax.tick_params(labelsize=16, bottom=True)
        ax.set_xlabel(
            'vector angle [degrees]', fontsize=16
        )
        ax.set_ylabel('rel. energy [kJ mol$^{-1}$]', fontsize=16)
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 21)

        fig.tight_layout()
        fig.savefig(
            f'{lig_name}_flex_torsion_energy.pdf',
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()


def plot_ladev_angle(json_files):
    # xbins = np.arange(0, 15, xwidth)
    fig, ax = plt.subplots(figsize=(8, 5))

    prop = {'c': 'gold', 'm': 'o'}

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
        ax.scatter(
            x=x,
            y=y,
            c=prop['c'],
            marker=prop['m'],
            alpha=1.0,
            edgecolor='k',
            s=160,
        )

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
    plot_ladev_planedev(json_files)
    plot_ladev_angle(json_files)
    plot_vectangle_energy_permol(json_files)


if __name__ == "__main__":
    main()
