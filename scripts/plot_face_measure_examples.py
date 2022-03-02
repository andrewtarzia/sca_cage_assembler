#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot example face measures.

Author: Andrew Tarzia

Date Created: 02 Mar 2022

"""

import os
import numpy as np
import matplotlib.pyplot as plt


def ii_model(a, c):
    aniso = np.linspace(1, 2, 100)
    b = aniso * a
    b
    mmdist1 = c+a
    mmdist2 = b

    x = []
    y = []
    for bb in mmdist2:
        relerror = (abs(bb-mmdist1)/max([mmdist1, bb]))*100
        x.append(bb/a)
        y.append(relerror)
    return x, y

def main():
    _figure_path = 'figures'

    c = 1
    _as = [1, 2, 3, 4]

    fig, ax = plt.subplots(figsize=(8, 5))
    for a in _as:
        x, y = ii_model(a, c)
        ax.plot(
            x, y,
            alpha=1.0,
            lw=2,
            label=f'a={a}'
        )

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('aspect ratio', fontsize=16)
    ax.set_ylabel('face mismatch [%]', fontsize=16)
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(
        os.path.join(_figure_path, 'face_examples.pdf'),
        dpi=720,
        bbox_inches='tight',
    )

    plt.close()


if __name__ == "__main__":
    main()
