#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot distributions of cage properties in set.

Author: Andrew Tarzia

Date Created: 12 Feb 2021

"""

import pandas as pd
import os
import sys
import matplotlib.pyplot as plt



def distribution_plot(df, col_name, sets_to_plot):
    _figure_path = 'figures'

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

    fig, ax = plt.subplots(figsize=(8, 5))

    _x_positions = 0
    _x_names = []
    forms = []
    does_not_form = []
    for cageset in sets_to_plot:
        print_name = sets_to_plot[cageset]
        set_df = df[df['cageset'] == cageset]
        _x_positions += 1
        _x_names.append((_x_positions, print_name))
        for i, row in set_df.iterrows():
            outcome = True if row['outcome'] == 1 else False
            y_val = row[col_name]
            if y_val is None:
                continue
            if outcome:
                forms.append((_x_positions, float(y_val)))
            else:
                does_not_form.append((_x_positions, float(y_val)))

    ax.scatter(
        x=[i[0] for i in does_not_form],
        y=[i[1] for i in does_not_form],
        c='gray',
        marker='o',
        alpha=1.0,
        s=40,
        label='does not form',
    )
    ax.scatter(
        x=[i[0] for i in forms],
        y=[i[1] for i in forms],
        c='#4691C3',
        edgecolors='k',
        marker='o',
        alpha=1.0,
        s=180,
        label='forms',
    )

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylabel(yprops[col_name][0], fontsize=16)
    ax.set_ylim(yprops[col_name][1])
    ax.set_xticks([i[0] for i in _x_names])
    ax.set_xticklabels([i[1] for i in _x_names])
    fig.legend(fontsize=16)
    fig.tight_layout()
    fig.savefig(
        os.path.join(_figure_path, f"distribution_{col_name}.pdf"),
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_set_energies(data, filename, sets_to_plot):

    fig, ax = plt.subplots(figsize=(8, 5))

    _x_positions = 0
    stab_energies = []
    _x_names = []
    for cageset in sets_to_plot:
        print_name = sets_to_plot[cageset]
        set_values = data[cageset]
        print(print_name, set_values)
        formed_symm = []
        all_energies = []
        for i in set_values:
            if set_values[i][1] is True:
                formed_symm.append((i, set_values[i]))
            all_energies.append(set_values[i][0])

        if cageset == 'cl1_quad2_1':
            other_energies = [
                (i-min(all_energies))*2625.5 for i in all_energies
            ]
            _x_positions += 1
            _x_names.append((_x_positions, print_name))
            stab_energies.append(-100)
            ax.scatter(
                x=[_x_positions for i in other_energies],
                y=other_energies,
                c='gray',
                s=40,
            )
        else:
            if len(formed_symm) != 1:
                raise ValueError(
                    'Missing something; there should be one True case.'
                )
            print(formed_symm)
            print(all_energies)
            stabilisation_energy = formed_symm[0][1][0] - min(all_energies)
            print(stabilisation_energy)
            stab_energies.append(stabilisation_energy*2625.5)
            other_energies = [
                (i-min(all_energies))*2625.5 for i in all_energies
            ]
            _x_positions += 1
            _x_names.append((_x_positions, print_name))
            ax.scatter(
                x=[_x_positions for i in other_energies],
                y=other_energies,
                c='gray',
                s=40,
            )

    ax.scatter(
        x=range(1, _x_positions+1),
        y=stab_energies,
        c='r',
        s=180,
    )

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('ligand', fontsize=16)
    ax.set_ylabel(r'rel. energy [kJmol$^{-1}$]', fontsize=16)
    # ax.set_xlim((0, 1))
    ax.set_ylim(-0.1, None)
    ax.set_xticks([i[0] for i in _x_names])
    ax.set_xticklabels([i[1] for i in _x_names])

    fig.tight_layout()
    fig.savefig(
        filename,
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def main():
    first_line = (
        'Usage: plot_set_energies.py'
    )
    if (not len(sys.argv) == 1):
        print(f"""
{first_line}
    """)
        sys.exit()
    else:
        pass

    _figure_path = 'figures'
    all_cage_properties = pd.read_csv('all_cage_csv_data.csv')

    # Define sets.
    sets_to_plot = {
        'cl1_quad2_1': 'quad2_1',
        'cl1_quad2_12': 'jd235',
        'cl1_quad2_8': 'jd257',
        'cl1_quad2_3': 'jd301',
        'cl1_quad2_16': 'jd326',
        'cl1_quad2_2': 'jd354',
        'cl1_quad2_5': 'jd370',
        'cl1_quad2_9': 'jd490',
    }
    target_cols = [
        'octop', 'rellsesum', 'minitors', 'lsesum',
        'maxcrplan', 'maxMLlength', 'porediam',
        'relformatione', 'maxfacemetalpd', 'maxintangledev',
        'm_cube_shape'
    ]
    for col_name in target_cols:
        distribution_plot(
            df=all_cage_properties,
            col_name=col_name,
            sets_to_plot=sets_to_plot,
        )

    # Get total energies from GFN.
    set_gfn_energies = {i: {} for i in sets_to_plot}
    for i, row in all_cage_properties.iterrows():
        setname = row['cageset']
        if setname not in sets_to_plot:
            continue
        print(setname)
        symm = row['symmetry']
        forms = True if row['outcome'] == 1 else False
        print(forms)
        ey_file = f'C_{setname}_{symm}_optc.ey'
        if not os.path.exists(ey_file):
            print(ey_file)
            continue
        with open(ey_file, 'r') as f:
            lines = f.readlines()
        total_energy = float(lines[0].rstrip())
        set_gfn_energies[setname][symm] = (total_energy, forms)

    # Get total energies from DFT.

    print(set_gfn_energies)
    # Plot relative energy cf. more stable symmetry of expt symmetry.
    plot_set_energies(
        data=set_gfn_energies,
        filename=os.path.join(_figure_path, 'set_energies_xtb.pdf'),
        sets_to_plot=sets_to_plot,
    )
    raise SystemExit('waiting on DFT')
    plot_set_energies(
        data=set_dft_energies,
        filename='set_energies_dft.pdf',
        sets_to_plot=sets_to_plot,
    )


if __name__ == "__main__":
    main()
