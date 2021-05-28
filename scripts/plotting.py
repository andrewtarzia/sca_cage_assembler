#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Plotting utilities module.

Author: Andrew Tarzia

Date Created: 15 Mar 2020

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors


def colors_i_like(palette=None):
    """
    A list of colours I like to choose from.

    palette options:
        None
        Base
        IBM
        Wong
        Tol
        CB_pairs

    """

    if palette is None:
        return [
            '#FA7268', '#F8A72A', '#DAF7A6', '#900C3F', '#6BADB0',
            '#DB869D', '#F6D973', 'mediumvioletred',
            'skyblue', 'gold', 'palegreen', 'coral',
        ]
    elif palette == 'Base':
        return [
            '#D81B60', '#1E88E5', '#FFC107', '#FE6100', '#004D40'
        ]
    elif palette == 'IBM':
        return [
            '#648FFF', '#785EF0', '#DC267F', '#FE6100', '#FFB000'
        ]
    elif palette == 'Wong':
        return [
            '#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442',
            '#0072B2', '#D55E00', '#CC79A7'
        ]
    elif palette == 'Tol':
        return [
            '#332288', '#117733', '#44AA99', '#88CCEE', '#DDCC77',
            '#CC6677', '#AA4499', '#882255',
        ]
    elif palette == 'CB_pairs':
        return [
            '#FFC20A', '#0C7BDC', '#994F00', '#006CD1', '#E1BE6A',
            '#40B0A6', '#E66100', '#5D3A9B', '#1AFF1A', '#4B0092',
            '#FEFE62', '#D35FB7', '#005AB5', '#DC3220', '#1A85FF',
            '#D41159',
        ]


def heatmap(
    data,
    row_labels,
    col_labels,
    ax=None,
    cbar_kw={},
    cbarlabel="",
    **kwargs
):
    """
    Create a heatmap from a numpy array and two lists of labels.

    From: https://matplotlib.org/3.1.1/gallery/
    images_contours_and_fields/image_annotated_heatmap.html
    #sphx-glr-gallery-images-contours-and-fields-image-
    annotated-heatmap-py

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is
        plotted.  If not provided, use current axes or create a new
        one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.
        Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.

    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

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

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)
    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    # ax.grid(which="minor", color="k", linestyle='-', linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(
    im,
    data=None,
    valfmt="{x:.2f}",
    textcolors=["black", "white"],
    threshold=None,
    na_points=None,
    **textkw
):
    """
    A function to annotate a heatmap.

    From: https://matplotlib.org/3.1.1/gallery/
    images_contours_and_fields/image_annotated_heatmap.html
    #sphx-glr-gallery-images-contours-and-fields-image-
    annotated-heatmap-py

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.
        Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should
        either use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used
        for values below a threshold, the second for those above.
        Optional.
    threshold
        Value in data units according to which the colors from
        textcolors are applied.  If None (the default) uses the middle
        of the colormap as separation.  Optional.
    na_points
        Tuple of (row, column) indices whose text should be manually
        set to N/A. Optional
    **kwargs
        All other arguments are forwarded to each call to `text` used
        to create the text labels.

    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(
        horizontalalignment="center",
        verticalalignment="center"
    )
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(
                color=textcolors[int(im.norm(data[i, j]) > threshold)]
            )
            cond1 = na_points is not None
            cond2 = False if na_points is None else (i, j) in na_points
            if cond1 and cond2:
                text = im.axes.text(j, i, 'N/A', **kw)
            else:
                text = im.axes.text(
                    j,
                    i,
                    valfmt(data[i, j], None),
                    **kw
                )
            texts.append(text)

    return texts


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0,
                    name='shiftedcmap'):
    """
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    From Stack Exchange:
        https://stackoverflow.com/questions/7404116/
        defining-the-midpoint-of-a-colormap-in-matplotlib

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
    """
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


def define_plot_cmap(fig, ax, mid_point, cmap, ticks, labels,
                     cmap_label):
    """
    Define cmap shifted to midpoint and plot colourbar

    """
    new_cmap = shiftedColorMap(
        cmap,
        midpoint=mid_point,
        name='shifted'
    )
    X = np.linspace(0, 1, 256)
    cax = ax.scatter(-X-100, -X-100, c=X, cmap=new_cmap)
    cbar = fig.colorbar(cax, ticks=ticks, spacing='proportional')
    cbar.ax.set_yticklabels(labels, fontsize=16)
    cbar.set_label(cmap_label, fontsize=16)
    return new_cmap


def parity_plot(X, Y, xtitle, ytitle, lim, c=None, marker=None):
    """
    Make parity plot.

    """
    if c is None:
        C = colors_i_like()[2]
    else:
        C = c
    if marker is None:
        M = 'o'
    else:
        M = marker
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(X, Y, c=C, edgecolors='k',
               marker=M, alpha=1.0, s=80)
    ax.plot(np.linspace(min(lim) - 1, max(lim) + 1, 2),
            np.linspace(min(lim) - 1, max(lim) + 1, 2),
            c='k', alpha=0.4)
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(xtitle, fontsize=16)
    ax.set_ylabel(ytitle, fontsize=16)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    return fig, ax


def histogram_plot_N(
    Y,
    X_range,
    width,
    alpha,
    color,
    edgecolor,
    xtitle,
    labels=None,
    density=False,
    N=1
):
    """
    Make histogram plot with 1 distribution.

    """

    fig, ax = plt.subplots(figsize=(8, 5))
    X_bins = np.arange(X_range[0], X_range[1], width)
    if N == 1:
        hist, bin_edges = np.histogram(
            a=Y,
            bins=X_bins,
            density=density
        )
        ax.bar(
            bin_edges[:-1],
            hist,
            align='edge',
            alpha=alpha,
            width=width,
            color=color,
            edgecolor=edgecolor
        )
    else:
        for i_ in range(N):
            if type(color) is not list or len(Y) != N:
                raise ValueError(
                    'Make sure color and Y are of length N'
                )
            hist, bin_edges = np.histogram(
                a=Y[i_],
                bins=X_bins,
                density=density
            )
            if labels[i_] is None:
                label = ''
            else:
                label = labels[i_]
            ax.bar(
                bin_edges[:-1],
                hist,
                align='edge',
                alpha=alpha[i_],
                width=width,
                color=color[i_],
                edgecolor=edgecolor,
                label=label
            )
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(xtitle, fontsize=16)
    if density is False:
        ax.set_ylabel('count', fontsize=16)
    elif density is True:
        ax.set_ylabel('frequency', fontsize=16)
    ax.set_xlim(X_range)
    if N > 1 and labels[0] is not None:
        ax.legend(fontsize=16)
    return fig, ax
