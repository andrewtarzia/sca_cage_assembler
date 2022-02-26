#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot categorisation of cage stability measures.

Author: Andrew Tarzia

Date Created: 12 Feb 2021

"""

import pandas as pd
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
from sklearn import tree


def main():
    first_line = (
        'Usage: decision_tree.py'
    )
    if (not len(sys.argv) == 1):
        print(f"""
{first_line}
    """)
        sys.exit()

    _figure_path = 'figures'

    all_cage_properties = pd.read_csv('all_cage_csv_data.csv')
    all_cage_properties = all_cage_properties.where(
        pd.notnull(all_cage_properties), None
    )

    target_cols = [
        # 'octop',
        # 'rellsesum',
        # 'minitors',
        # 'maxcrplan',
        'maxMLlength',
        # 'relformatione',
        # 'maxfacemetalpd',
        'maxintangledev',
        'm_cube_shape'
    ]

    X = []
    y = []
    for i, row in all_cage_properties.iterrows():
        if row['m_cube_shape'] is None:
            continue
        x = []
        flagged = False
        for tc in target_cols:
            try:
                x.append(float(row[tc ]))
            except TypeError:
                flagged = True

        if not flagged:
            X.append(x)
            y.append(int(row['outcome']))

    X = np.array(X)
    print(X.shape)
    y = np.array(y)
    print(y.shape)

    clf = tree.DecisionTreeClassifier()
    clf = clf.fit(X, y)

    fig, ax = plt.subplots(figsize=(8, 8))
    tree.plot_tree(
        decision_tree=clf,
        ax=ax,
        feature_names=target_cols,
        class_names=['no', 'yes']
    )

    fig.tight_layout()
    fig.savefig(
        os.path.join(_figure_path, 'decision_tree.pdf'),
        dpi=160,
        bbox_inches='tight'
    )
    plt.close()



if __name__ == "__main__":
    main()
