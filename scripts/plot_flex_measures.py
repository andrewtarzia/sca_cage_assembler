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


def plot_ladists(json_files):
    num_mols = len(json_files)
    xwidth = 0.2
    xbins = np.arange(7.5, 25, xwidth)
    # xbins = np.arange(0, 15, xwidth)
    fig, ax = plt.subplots(num_mols, 1, sharex=True, figsize=(8, 10))
    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=0)

    for i, _file in enumerate(json_files):
        lig_name = (
            _file.replace('_flex_measure.json', '').split('_')[1]
        )
        with open(_file, 'r') as f:
            data = json.load(f)
        measures = data['long_axis_distances']
        # measures = [i-min(measures) for i in measures]
        ax[i].text(7.5, 0.5, lig_name, fontsize=16)
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
    ax[i].set_xlabel(
        r'long-axis deviation [$\mathrm{\AA}$]', fontsize=16
    )
    ax[i].set_ylabel('frequency', fontsize=16)
    fig.tight_layout()
    fig.savefig(
        'flex_dists.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_ladev_planedev(json_files):
    # xbins = np.arange(0, 15, xwidth)
    fig, ax = plt.subplots(figsize=(8, 5))
    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=0)

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
        ax.scatter(
            x=x,
            y=y,
            c='gold',
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


def main():
    first_line = (
        'Usage: plot_flex_measures.py'
    )
    if (not len(sys.argv) == 1):
        print(f"""{first_line}""")
        sys.exit()

    json_files = glob('*_flex_measure.json')

    plot_ladists(json_files)
    plot_ladev_planedev(json_files)


if __name__ == "__main__":
    main()
