#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot example shape measures.

Author: Andrew Tarzia

Date Created: 08 Nov 2020

"""

from glob import glob
import matplotlib.pyplot as plt
import stk

from utilities import calculate_cube_shape_measure


def main():

    structure_files = sorted(glob('*.mol'))
    print(structure_files)

    shapes = {}
    for fi in structure_files:
        print(fi)
        name = fi.replace('.mol', '')
        struct = stk.BuildingBlock.init_from_file(fi)
        print(struct)
        shapes[name] = calculate_cube_shape_measure(
            name=name,
            molecule=struct,
        )['CU-8']

    print(shapes)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(
        [int(i.replace('shape_', '')) for i in shapes],
        [shapes[i] for i in shapes],
        c='gold',
        edgecolors='k',
        marker='o',
        alpha=1.0,
        s=160,
    )

    # Set number of ticks for x-axis
    x_pos = [int(i.replace('shape_', '')) for i in shapes]
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xticks(x_pos)
    ax.set_xticklabels([str(i) for i in x_pos])
    ax.set_xlim(0, max(x_pos)+1)
    ax.set_xlabel('structure', fontsize=16)
    ax.set_ylabel('CU-8 shape measure', fontsize=16)

    fig.tight_layout()
    fig.savefig('shape_examples.pdf', dpi=720, bbox_inches='tight')

    plt.close()


if __name__ == "__main__":
    main()
