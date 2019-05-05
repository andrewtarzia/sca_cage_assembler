#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to show the effect of the LHS and RHS changes on the trapezoid geometry

Author: Andrew Tarzia

Date Created: 05 May 2019

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
sys.path.insert(0, '/home/atarzia/thesource/')
from calculations import angle_between


def plot_trapezoid(patches, pts, c):
    '''Plot trapezoids.

    '''
    p = Polygon(pts, closed=True)
    patches.append(p)
    return patches


def main():
    '''Build trapzeoids to plot.

    '''
    fig, ax = plt.subplots(figsize=(5, 5))
    # create patches and add LHS, RHS on figure
    patches = []
    PTS = np.array([[5, 4], [5, 16],
                    [15, 15], [15, 5]])
    NN1 = PTS[1][1] - PTS[0][1]
    NN2 = PTS[2][1] - PTS[3][1]
    print(NN1, NN2)
    LHS = round((NN1 - NN2) / 2, 2)
    v1 = PTS[2] - PTS[1]
    v2 = PTS[2] - PTS[3]
    print(v1, v2)
    angle = np.degrees(angle_between(v1, v2))
    print(angle)
    L = np.linalg.norm(v1)
    print(L)
    RHS = round(L * np.cos(np.radians(180 - angle)), 2)
    print(RHS)
    ax.text(1, 18, str(LHS)+' = '+str(RHS))
    v1 = PTS[3] - PTS[2]
    v2 = PTS[3] - PTS[0]
    print(v1, v2)
    angle = np.degrees(angle_between(v1, v2))
    print(angle)
    L = np.linalg.norm(v2)
    print(L/2)
    RHS = round(L * np.cos(np.radians(180 - angle)), 2)
    print(RHS)
    ax.text(1, 1, str(LHS)+' = '+str(RHS))
    patches = plot_trapezoid(patches=patches, pts=PTS, c='k')

    # set up frame
    p = PatchCollection(patches, alpha=1.0, color=None, edgecolor='k')
    ax.add_collection(p)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('X', fontsize=16)
    ax.set_ylabel('Y', fontsize=16)
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 20)
    fig.tight_layout()
    fig.savefig('viz_trapezoid.pdf', dpi=720,
                bbox_inches='tight')


if __name__ == "__main__":
    main()
