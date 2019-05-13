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
import calculations


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
    a_line_pt = [5, 5]
    b_line_pt = [15, 5]
    # set vectors and angles
    A = 10
    B = 9
    C = 2.06
    ABon2 = A-B / 2
    cospimind = 2*(C)/ABon2
    cospiminc = 2*(C)/ABon2
    d = np.degrees(np.arccos(cospimind) + np.pi)
    c = np.degrees(np.arccos(cospiminc) + np.pi)
    a = 180 - d
    b = 180 - c
    # define trapezoid points based on variables
    pt1 = [a_line_pt[0], a_line_pt[1]+A/2]
    pt2 = [b_line_pt[0], b_line_pt[1]+B/2]
    pt3 = [b_line_pt[0], b_line_pt[1]-B/2]
    pt4 = [a_line_pt[0], a_line_pt[1]-A/2]
    PTS = np.array([pt1, pt2, pt3, pt4])
    NN1 = PTS[1][1] - PTS[0][1]
    NN2 = PTS[2][1] - PTS[3][1]
    print(NN1, NN2)
    LHS = round((NN1 - NN2) / 2, 2)
    v1 = PTS[2] - PTS[1]
    v2 = PTS[2] - PTS[3]
    print(v1, v2)
    angle = np.degrees(calculations.angle_between(v1, v2))
    print(angle)
    L = np.linalg.norm(v1)
    print(L)
    RHS = round(L * np.cos(np.radians(180 - angle)), 2)
    print(RHS)
    v1 = PTS[3] - PTS[2]
    v2 = PTS[3] - PTS[0]
    print(v1, v2)
    angle = np.degrees(calculations.angle_between(v1, v2))
    print(angle)
    L = np.linalg.norm(v2)
    print(L/2)
    RHS = round(L * np.cos(np.radians(180 - angle)), 2)
    print(RHS)
    ax.text(1, 1, 'A ='+str(A)+' B ='+str(B))
    ax.text(1, 1, 'a ='+str(a)+' b ='+str(b)+' c ='+str(c)+' d ='+str(d))
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
