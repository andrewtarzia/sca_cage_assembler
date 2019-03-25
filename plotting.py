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
    ax.plot(np.linspace(-1, max(lim)+1, 2),
            np.linspace(-1, max(lim)+1, 2),
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
