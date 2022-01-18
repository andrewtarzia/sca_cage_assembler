#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot categorisation of cage stability measures.

Author: Andrew Tarzia

Date Created: 12 Feb 2021

"""

import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
from sklearn.datasets import load_iris
from sklearn import tree


def categorisation_plot(df, xray_df, col_name):

    yprops = {
        'octop': (
            r'min. $q_{\mathrm{oct}}$', (None, None), 'min'
        ),
        'm_cube_shape': ('CU-8 cube measure', (-0.5, None), 'min'),
        'rellsesum': (
            r'rel. sum strain energy [kJmol$^{-1}$]',
            (-10, 1000),
            'max'
        ),
        'lsesum': (
            r'sum strain energy [kJmol$^{-1}$]',
            (None, None),
            'max'
        ),
        'minitors': (
            r'min. imine torsion [degrees]', (None, None), 'min'
        ),
        'maxcrplan': (
            r'max. ligand distortion [$\mathrm{\AA}$]',
            (None, None),
            'max'
        ),
        'maxdifffaceaniso': (
            r'max. $\Delta$opposing face anisotropy [%]',
            (None, None),
            'max'
        ),
        'maxMLlength': (
            r'max. N-Zn bond length [$\mathrm{\AA}$]',
            (None, None),
            'max'
        ),
        'porediam': (
            r'pore diameter [$\mathrm{\AA}$]', (None, None), 'min'
        ),
        'relformatione': (
            r'rel. formation energy [kJmol$^{-1}$]', (-10, 1000), 'max'
        ),
        'maxfacemetalpd': (
            r'max. face metal planarity deviation [$\mathrm{\AA}$]',
            (None, None),
            'max'
        ),
        'maxintangledev': (
            r'max. interior angle deviation [degrees]',
            (None, None),
            'max'
        ),
    }

    dont_show = ['rellsesum', 'relformatione']

    dx = 0.15
    fig, ax = plt.subplots(figsize=(5, 5))

    forms_df = df[df['outcome'] == 1]
    does_not_form_df = df[df['outcome'] == 0]

    forms = list(forms_df[col_name])
    does_not_form = list(does_not_form_df[col_name])
    xray = list(xray_df[col_name])
    for i in forms:
        ax.scatter(
            0.25+(dx*(np.random.random() - 0.5) * 2),
            i,
            c='#4691C3',
            edgecolors='k',
            marker='o',
            alpha=1.0,
            s=120,
        )

    if col_name not in dont_show:
        for i in xray:
            ax.scatter(
                0.25+(dx*(np.random.random() - 0.5) * 2),
                i,
                c='#E74C3C',
                edgecolors='k',
                marker='X',
                alpha=1.0,
                s=120,
            )

    for i in does_not_form:
        ax.scatter(
            0.75+(dx*(np.random.random() - 0.5) * 2),
            i,
            c='#B8BEC3',
            edgecolors='k',
            marker='o',
            alpha=1.0,
            s=60,
        )

    if yprops[col_name][2] == 'max':
        if max(forms) < min(does_not_form):
            ax.axhline(y=max(forms), c='k', lw=2, linestyle='--')
    elif yprops[col_name][2] == 'min':
        if min(forms) > max(does_not_form):
            ax.axhline(y=min(forms), c='k', lw=2, linestyle='--')

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylabel(yprops[col_name][0], fontsize=16)
    ax.set_xlim((0, 1))
    ax.set_ylim(yprops[col_name][1])
    ax.set_xticks([0.25, 0.75])
    ax.set_xticklabels(['forms', 'does not form'])

    if col_name not in dont_show:
        ax.scatter(
            -100, 100,
            c='none',
            edgecolors='k',
            marker='o',
            alpha=1.0,
            s=120,
            label='calculated structure',
        )
        ax.scatter(
            -100, 100,
            c='#E74C3C',
            edgecolors='k',
            marker='X',
            alpha=1.0,
            s=120,
            label='xray structure',
        )
        ax.legend(fontsize=16)
    fig.tight_layout()
    fig.savefig(
        f"categorical_{col_name}.pdf",
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def main():
    first_line = (
        'Usage: decision_tree.py'
    )
    if (not len(sys.argv) == 1):
        print(f"""
{first_line}
    """)
        sys.exit()

    all_cage_properties = pd.read_csv('all_cage_csv_data.csv')
    all_cage_properties = all_cage_properties.where(
        pd.notnull(all_cage_properties), None
    )

    target_cols = [
        # 'octop',
        # 'rellsesum',
        # 'minitors',
        # 'maxcrplan',
        # 'maxdifffaceaniso',
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
        'decision_tree.pdf',
        dpi=160,
        bbox_inches='tight'
    )
    plt.close()



if __name__ == "__main__":
    main()
