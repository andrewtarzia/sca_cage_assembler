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
    a_line_pt = [6, 10]
    b_line_pt = [10, 10]
    # set vectors and angles
    A = 10
    B = 7
    # define trapezoid points based on variables
    pt1 = np.asarray([a_line_pt[0], a_line_pt[1]+A/2])
    pt2 = np.asarray([b_line_pt[0], b_line_pt[1]+B/2])
    pt3 = np.asarray([b_line_pt[0], b_line_pt[1]-B/2])
    pt4 = np.asarray([a_line_pt[0], a_line_pt[1]-A/2])
    C = np.linalg.norm(pt1-pt2) / 2
    print('C', C)
    ABon2 = (A-B) / 2
    print('ABon2', ABon2)
    cospimind = ABon2/(2*(C))
    cospiminc = ABon2/(2*(C))
    print('cos', cospiminc, cospimind)
    print('acos', np.arccos(cospiminc), np.arccos(cospimind))
    d = np.degrees(np.pi - np.arccos(cospimind))
    c = np.degrees(np.pi - np.arccos(cospiminc))
    a = 180 - d
    b = 180 - c
    print('angles:', a, b, c, d)
    PTS = np.array([pt1, pt2, pt3, pt4])
    NN1 = A
    NN2 = B
    LHS = round((NN1 - NN2) / 2, 2)
    v1 = PTS[1] - PTS[0]
    v2 = PTS[1] - PTS[2]
    angle = np.degrees(calculations.angle_between(v1, v2))
    print('d:', angle)
    L = np.linalg.norm(v1)
    print('2C:', L)
    RHS = round(L * np.cos(np.radians(180 - angle)), 2)
    print('RHS', RHS)
    ax.text(1, 15, 'A = '+str(round(A, 2)))
    ax.text(1, 13, 'B = '+str(round(B, 2)))
    ax.text(1, 11, 'C = '+str(round(C, 2)))
    ax.text(1, 9, 'a = '+str(round(a, 2)))
    ax.text(1, 7, 'b = '+str(round(b, 2)))
    ax.text(1, 5, 'c = '+str(round(c, 2)))
    ax.text(1, 3, 'd = '+str(round(d, 2)))
    ax.text(1, 1, str(round(LHS, 2))+' = '+str(round(RHS, 2)))
    patches = plot_trapezoid(patches=patches, pts=PTS, c='k')

    # set up frame
    p = PatchCollection(patches, alpha=1.0, color=None, edgecolor='k')
    ax.add_collection(p)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('X', fontsize=16)
    ax.set_ylabel('Y', fontsize=16)
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 20)
    fig.tight_layout()
    fig.savefig('viz_trapezoid.pdf', dpi=720,
                bbox_inches='tight')


if __name__ == "__main__":
    main()
