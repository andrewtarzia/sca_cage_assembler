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
import matplotlib.colors as colors


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
    plt.close()


def scatter_plot(X, Y, outfile, xtitle, ytitle, xlim, ylim, title=None,
                 c='firebrick', edgecolors='k',
                 marker='o', alpha=1.0, s=80, Z=None, cmap=None):
    '''Make scatter plot.

    '''
    fig, ax = plt.subplots(figsize=(8, 5))
    if cmap is None and Z is None:
        ax.scatter(X, Y, c=c, edgecolors=edgecolors,
                   marker=marker, alpha=alpha, s=s)
    else:
        cmp = define_plot_cmap(fig, ax,
                               mid_point=cmap['mid_point'],
                               cmap=cmap['cmap'],
                               ticks=cmap['ticks'],
                               labels=cmap['labels'],
                               cmap_label=cmap['cmap_label'])
        ax.scatter(X, Y, c=cmp(Z), edgecolors=edgecolors,
                   marker=marker, alpha=alpha, s=s)
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(xtitle, fontsize=16)
    ax.set_ylabel(ytitle, fontsize=16)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    if title is not None:
        ax.set_title(title, fontsize=16)
    fig.tight_layout()
    fig.savefig(outfile, dpi=720,
                bbox_inches='tight')
    plt.close()


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
    plt.close()


def flat_line(ax, x, y, w=0, C='k', m='x'):
    ax.plot([x - w, x, x + w], [y, y, y], c=C)
    ax.scatter(x, y, marker=m, c=C)


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    From Stack Exchange:
        https://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }
    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)
    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)
        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))
    newcmap = colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)
    return newcmap


def define_plot_cmap(fig, ax, mid_point, cmap, ticks, labels, cmap_label):
    """Define cmap shifted to midpoint and plot colourbar

    """
    new_cmap = shiftedColorMap(cmap, midpoint=mid_point, name='shifted')
    X = np.linspace(0, 1, 256)
    cax = ax.scatter(-X-100, -X-100, c=X, cmap=new_cmap)
    cbar = fig.colorbar(cax, ticks=ticks, spacing='proportional')
    cbar.ax.set_yticklabels(labels, fontsize=16)
    cbar.set_label(cmap_label, fontsize=16)
    return new_cmap
