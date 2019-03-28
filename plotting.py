#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for base plots.

Author: Andrew Tarzia

Date Created: 25 Mar 2019
"""

import numpy as np
import matplotlib.pyplot as plt


def parity_plot(X, Y, outfile, xtitle, ytitle, lim):
    '''Make parity plot.

    '''
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(X, Y, c='firebrick', edgecolors='k',
               marker='o', alpha=1.0, s=80)
    ax.plot(np.linspace(min(lim) - 1, max(lim) + 1, 2),
            np.linspace(min(lim) - 1, max(lim) + 1, 2),
            c='k', alpha=0.4)
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(xtitle, fontsize=16)
    ax.set_ylabel(ytitle, fontsize=16)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    fig.tight_layout()
    fig.savefig(outfile, dpi=720,
                bbox_inches='tight')


def histogram_plot_1(Y, X_range, width, alpha, color, edgecolor,
                     outfile, xtitle, density=False):
    '''Make histogram plot with 1 distribution.

    '''

    fig, ax = plt.subplots(figsize=(8, 5))
    X_bins = np.arange(X_range[0], X_range[1], width)
    hist, bin_edges = np.histogram(a=Y, bins=X_bins, density=density)
    ax.bar(bin_edges[:-1],
           hist,
           align='edge',
           alpha=alpha, width=width,
           color=color,
           edgecolor=edgecolor)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(xtitle, fontsize=16)
    if density is False:
        ax.set_ylabel('count', fontsize=16)
    elif density is True:
        ax.set_ylabel('frequency', fontsize=16)
    ax.set_xlim(X_range)
    fig.tight_layout()
    fig.savefig(outfile, dpi=720,
                bbox_inches='tight')


def flat_line(ax, x, y, w=0, C='k', m='x'):
    ax.plot([x - w, x, x + w], [y, y, y], c=C)
    ax.scatter(x, y, marker=m, c=C)
