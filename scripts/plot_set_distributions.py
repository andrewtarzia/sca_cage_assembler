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
        cset_ys = []
        for i, row in set_df.iterrows():
            outcome = True if row['outcome'] == 1 else False
            y_val = row[col_name]
            if y_val is None:
                continue
            if outcome:
                forms.append((_x_positions, float(y_val)))
            else:
                does_not_form.append((_x_positions, float(y_val)))
            cset_ys.append(float(y_val))

        parts = ax.violinplot(
            cset_ys,
            [_x_positions],
            # points=200,
            vert=True,
            widths=0.8,
            showmeans=False,
            showextrema=False,
            showmedians=False,
            bw_method=0.5,
        )

        for pc in parts['bodies']:
            pc.set_facecolor('gray')
            pc.set_edgecolor('none')
            pc.set_alpha(0.3)


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
        c='gold',
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
        formed_symm = []
        all_energies = []
        for i in set_values:
            if set_values[i][1] is True:
                formed_symm.append((i, set_values[i][0]))
            all_energies.append(set_values[i][0])

        if cageset == 'cl1_quad2_1':
            other_energies = [
                (i-min(all_energies))*2625.5 for i in all_energies
            ]
            _x_positions += 1
            _x_names.append((_x_positions, print_name))
            stab_energies.append(-100)

        else:
            if len(formed_symm) < 1:
                raise ValueError(
                    'Missing something; there should be one True case.'
                )

            other_energies = [
                (i-min(all_energies))*2625.5 for i in all_energies
            ]
            _x_positions += 1
            _x_names.append((_x_positions, print_name))
            for i in formed_symm:
                print(i)
                stabilisation_energy = (
                    i[1] - min(all_energies)
                )
                stab_energies.append(
                    (_x_positions, stabilisation_energy*2625.5)
                )
                print(stab_energies)

        parts = ax.violinplot(
            other_energies,
            [_x_positions],
            # points=200,
            vert=True,
            widths=0.8,
            showmeans=False,
            showextrema=False,
            showmedians=False,
            bw_method=0.5,
        )

        for pc in parts['bodies']:
            pc.set_facecolor('gray')
            pc.set_edgecolor('none')
            pc.set_alpha(0.3)
        ax.scatter(
            x=[_x_positions for i in other_energies],
            y=other_energies,
            c='gray',
            s=40,
            alpha=0.5,
        )

    ax.scatter(
        x=[i[0] for i in stab_energies],
        y=[i[1] for i in stab_energies],
        c='gold',
        edgecolors='k',
        s=180,
    )

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('tetra-aniline', fontsize=16)
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


def plot_set_parities(dftdata, xtbdata, filename, sets_to_plot):

    fig, ax = plt.subplots(figsize=(5, 5))

    for cageset in sets_to_plot:
        if cageset != 'cl1_quad2_12':
            continue
        xtb_set_values = xtbdata[cageset]
        dft_set_values = dftdata[cageset]

        all_energies = []
        for i in xtb_set_values:
            if xtb_set_values[i][1] is True:
                formed_pair = (
                    (xtb_set_values[i][0], dft_set_values[i][0])
                )
            all_energies.append(
                (xtb_set_values[i][0], dft_set_values[i][0])
            )

        minx_energy = min([i[0] for i in all_energies])
        miny_energy = min([i[1] for i in all_energies])
        ax.scatter(
            x=[
                (i[0]-minx_energy)*2625.5 for i in all_energies
            ],
            y=[
                (i[1]-miny_energy)*2625.5 for i in all_energies
            ],
            c='gray',
            s=40,
            alpha=0.5,
        )
        ax.scatter(
            x=[(formed_pair[0]-minx_energy)*2625.5],
            y=[(formed_pair[1]-miny_energy)*2625.5],
            c='gold',
            edgecolors='k',
            s=180,
        )

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(r'rel. GFN2-xTB energy [kJmol$^{-1}$]', fontsize=16)
    ax.set_ylabel(r'rel. DFT energy [kJmol$^{-1}$]', fontsize=16)
    ax.set_xlim(0, 600)
    ax.set_ylim(0, 600)

    fig.tight_layout()
    fig.savefig(
        filename,
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def main():
    first_line = (
        'Usage: plot_set_distributions.py'
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
    all_cage_properties = all_cage_properties.where(
        pd.notnull(all_cage_properties), None
    )

    # Define sets.
    sets_to_plot = {
        'cl1_quad2_5': 'A',
        'cl1_quad2_16': 'B',
        'cl1_quad2_12': 'C',
        'cl1_quad2_3': 'D',
        'cl1_quad2_8': 'E',
        'cl1_quad2_2': 'F',
        # 'cl1_quad2_1': 'G',
    }
    target_cols = [
        'octop', 'rellsesum', 'minitors', 'lsesum',
        'maxcrplan', 'maxMLlength', 'porediam',
        'relformatione', 'maxintangledev',
        'm_cube_shape'
    ]
    for col_name in target_cols:
        distribution_plot(
            df=all_cage_properties,
            col_name=col_name,
            sets_to_plot=sets_to_plot,
        )


    set_gfn_energies = {i: {} for i in sets_to_plot}
    set_dft_energies = {i: {} for i in sets_to_plot}
    for i, row in all_cage_properties.iterrows():
        setname = row['cageset']
        if setname not in sets_to_plot:
            continue
        symm = row['symmetry']
        forms = True if row['outcome'] == 1 else False

        # Get total energies from GFN.
        ey_file = f'C_{setname}_{symm}_optc.ey'
        if not os.path.exists(ey_file):
            print(ey_file, 'missing')
            continue
        with open(ey_file, 'r') as f:
            lines = f.readlines()
        total_energy = float(lines[0].rstrip())
        set_gfn_energies[setname][symm] = (total_energy, forms)

        # Get total energies from DFT.
        ey_file = os.path.join(
            'set_dft_run', f'C_{setname}_{symm}_optc_opt.out'
        )
        if not os.path.exists(ey_file):
            print(ey_file, 'missing')
            continue
        with open(ey_file, 'r') as f:
            for line in f.readlines():
                if (
                    'ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:'
                ) in line:
                    energy_line = line
                    break

        total_energy = float(energy_line.strip().split()[-1])
        set_dft_energies[setname][symm] = (total_energy, forms)

    # Plot relative energy cf. more stable symmetry of expt symmetry.
    plot_set_energies(
        data=set_gfn_energies,
        filename=os.path.join(_figure_path, 'set_energies_xtb.pdf'),
        sets_to_plot=sets_to_plot,
    )
    print(set_gfn_energies)
    print('---')
    plot_set_parities(
        dftdata=set_dft_energies,
        xtbdata=set_gfn_energies,
        filename=os.path.join(_figure_path, 'set_energies_dft.pdf'),
        sets_to_plot=sets_to_plot,
    )
    print(set_dft_energies)


if __name__ == "__main__":
    main()
