#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build complex library.

Author: Andrew Tarzia

Date Created: 27 Jan 2020
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

import cage_building
from utilities import read_lib


def plot_set_Y_vs_aniso(data, Y_name, ylabel, ylim, filename):
    C = '#AFE074'
    M = 'o'

    fig, ax = plt.subplots(figsize=(8, 5))
    for name in data:
        X = data[name]['aspect_ratio']
        # Iterate over all cages in set.
        for Y_d in data[name][Y_name]:
            Y = data[name][Y_name][Y_d]
            print(name, X, Y)
            print(data[name])
            print(data[name][Y_name])
            input()
            ax.scatter(
                X,
                Y,
                c=C,
                edgecolors='k',
                marker=M,
                alpha=1.0,
                s=120
            )
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('ligand aspect ratio [1:X]', fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_xlim(1, 10)
    ax.set_ylim(ylim)

    fig.tight_layout()
    fig.savefig(
        filename,
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def het_prism_analysis(cage_set):
    """
    Analyse cage set.

    Analysis performed:
         -
         -

    """

    cage_set.load_properties()
    # Compare average pore volume of each of the three topologies.
    three_top = {'m4l4spacer': [], 'm8l6face': [], 'm6l2l3': []}
    built_prop = cage_set.built_cage_properties
    all_pore_volumes = []
    all_min_OPs = []
    all_topo_strs = []
    for C in cage_set.cages_to_build:
        C_data = built_prop[C.name]
        TOPO = C.topology_string
        three_top[TOPO].append(
            C_data['pw_prop']['pore_volume_opt']
        )
        all_pore_volumes.append(
            C_data['pw_prop']['pore_volume_opt']
        )
        all_topo_strs.append(TOPO)
    avg_pore_vol = {
        i: np.mean(three_top[i]) for i in three_top
    }
    print('avg pore volumes:', avg_pore_vol)

    # Ensure at least one prismatic cage is stable.
    prism_oct_op = {}
    for C in cage_set.cages_to_build:
        C_data = built_prop[C.name]
        TOPO = C.topology_string
        # Get minimium octahedral OP of the metal that is in the
        # complex building block.
        print(cage_set.complex_dicts)
        atom_no_of_interest = list(set([
            int(cage_set.complex_dicts[i]['metal_atom_no'])
            for i in cage_set.complex_dicts
        ]))
        print(atom_no_of_interest)
        print('rrrr',)
        C_OP = [
            C_data['op_prop'][str(i)]
            for i in atom_no_of_interest
        ]
        print(C_OP)
        target_OPs = [
            i[j]['oct']
            for i in C_OP
            for j in i
        ]
        print('t', target_OPs)
        all_min_OPs.append(min(target_OPs))
        if TOPO == 'm6l2l3':
            prism_oct_op[C.name] = min(target_OPs)
    print('minimum prism OPs:', prism_oct_op)

    # Plot all order parameter minimums VS average pore volumes.
    cage_set.plot_min_OPs_avg_PV(
        X=all_pore_volumes,
        Y=all_min_OPs,
        T=all_topo_strs
    )
    sys.exit()


def homo_cube_analysis(cage_set):
    """
    Analyse cage set.

    Analysis performed:
         -
         -

    """

    cage_set.load_properties()
    built_prop = cage_set.built_cage_properties

    # Get measures of all cages.
    measures = {
        'oct_op': {},
        'lse_max': {},
        'min_imine_torsions': {},
        'max_ligand_distortion': {},
        'max_diff_face_aniso': {},
    }
    for C in cage_set.cages_to_build:
        C_data = built_prop[C.name]
        # Get minimium octahedral OP of the metal that is in the
        # complex building block.
        atom_no_of_interest = list(set([
            int(cage_set.complex_dicts[i]['metal_atom_no'])
            for i in cage_set.complex_dicts
        ]))
        C_OP = [
            C_data['op_prop'][str(i)]
            for i in atom_no_of_interest
        ]
        target_OPs = [
            i[j]['oct']
            for i in C_OP
            for j in i
        ]
        measures['oct_op'][C.name] = min(target_OPs)
        measures['lse_max'][C.name] = max([
            C_data['li_prop']['strain_energies'][i]
            for i in C_data['li_prop']['strain_energies']
        ])
        measures['min_imine_torsions'][C.name] = min([
            j
            for i in C_data['li_prop']['imine_torsions']
            for j in C_data['li_prop']['imine_torsions'][i]
        ])
        measures['max_ligand_distortion'][C.name] = max([
            C_data['li_prop']['core_planarities'][i]
            for i in C_data['li_prop']['core_planarities']
        ])
        measures['max_diff_face_aniso'][C.name] = max([
            (i[2] - i[3]) / i[2]
            for i in C_data['fa_prop']
        ])

    print('minimum OPs:', measures['oct_op'])
    print('max LSE:', measures['lse_max'])
    print('min imine torsion:', measures['min_imine_torsions'])
    print('max core planarities:', measures['max_ligand_distortion'])
    print('max diff in face aniso:', measures['max_diff_face_aniso'])
    input('----------------------------------------------------')
    # Plot all order parameter minimums.
    cage_set.plot_Y(
        data=measures['oct_op'],
        ylabel=r'min. $q_{\mathrm{oct}}$',
        ylim=(0, 1),
        filename=f'{cage_set.name}_minOPs.pdf'
    )
    cage_set.plot_Y(
        data={
            i: measures['lse_max'][i]-min(measures['lse_max'].values())
            for i in measures['lse_max']
        },
        ylabel=r'rel. max. strain energy [kJ/mol]',
        ylim=(-4, 50),
        filename=f'{cage_set.name}_maxLSE.pdf'
    )
    cage_set.plot_Y(
        data=measures['min_imine_torsions'],
        ylabel=r'min. imine torsion [degrees]',
        ylim=(0, 185),
        filename=f'{cage_set.name}_mintors.pdf'
    )
    cage_set.plot_Y(
        data=measures['max_ligand_distortion'],
        ylabel=r'max. ligand distortion [$\mathrm{\AA}$]',
        ylim=(0, 185),
        filename=f'{cage_set.name}_maxdistortion.pdf'
    )
    cage_set.plot_Y(
        data=measures['max_diff_face_aniso'],
        ylabel=(
            r'max. $\Delta$opposing face anisotropy [%]'
        ),
        ylim=(-0.1, 0.5),
        filename=f'{cage_set.name}_maxfadiff.pdf'
    )

    return measures


def analyse_cages(cage_sets):

    AR_data = {}
    for cage_set in cage_sets:
        if isinstance(cage_set, cage_building.HetPrism):
            het_prism_analysis(cage_set)
        if isinstance(cage_set, cage_building.HoCube):
            measures = homo_cube_analysis(cage_set)
            # For all HoCube sets, collect ligand aspect ratio data.
            AR_data[cage_set.name] = {
                'min_OPs': measures['oct_op'],
                'lse_max': {
                    i: (
                        measures['lse_max'][i] -
                        min(measures['lse_max'].values())
                    )
                    for i in measures['lse_max']
                },
                'min_tor': measures['min_imine_torsions'],
                'max_dis': measures['max_ligand_distortion'],
                'aspect_ratio': cage_set.ligand_aspect_ratio,
                'max_rfa': measures['max_diff_face_aniso'],
            }

    # Plot ligand aspect ratio data.
    tests = {
        # Test: (ylabel, ylim)
        'min_OPs': (r'min. $q_{\mathrm{oct}}$', (0, 1)),
        'lse_max': (r'rel. max. strain energy [kJ/mol]', (-4, 50)),
        'min_tor': (r'min. imine torsion [degrees]', (0, 185)),
        'max_dis': (
            r'max. ligand distortion [$\mathrm{\AA}$]', (0, 185)
        ),
        'max_rfa': (
            r'max. $\Delta$opposing face anisotropy [%]', (-0.1, 0.5)
        ),
    }
    if len(AR_data) > 0:
        for t in tests:
            plot_set_Y_vs_aniso(
                data=AR_data,
                Y_name=t,
                ylabel=tests[t][0],
                ylim=tests[t][0],
                filename=f'plotset_{t}_VA.pdf'
            )


def build_cages(
    ligands,
    complexes,
    cage_lib,
    ligand_directory,
    complex_directory
):

    cage_sets = []
    for name in cage_lib:
        cage_c = cage_lib[name]
        compl_names = cage_c['corners']
        comps = {i: complexes[i] for i in compl_names}

        if cage_c['heteroleptic']:
            cage_set = cage_building.HetPrism(
                name=name,
                cage_dict=cage_c,
                complex_dicts=comps,
                ligand_dicts=ligands,
                ligand_dir=ligand_directory,
                complex_dir=complex_directory
            )
        else:
            cage_set = cage_building.HoCube(
                name=name,
                cage_dict=cage_c,
                complex_dicts=comps,
                ligand_dicts=ligands,
                ligand_dir=ligand_directory,
                complex_dir=complex_directory
            )

        for C in cage_set.cages_to_build:
            C.build()

            default_free_e = C.free_electron_options[0]
            # Use a slightly different collapser threshold for
            # different topologies.
            if C.topology_string == 'm6l2l3':
                step_size = 0.05
                distance_cut = 3.0
                scale_steps = True
                expected_ligands = 2
            elif C.topology_string == 'm8l6face':
                step_size = 0.05
                distance_cut = 2.5
                scale_steps = False
                expected_ligands = 1
            else:
                step_size = 0.05
                distance_cut = 2.0
                scale_steps = True
                expected_ligands = 1

            C.optimize(
                free_e=default_free_e,
                step_size=step_size,
                distance_cut=distance_cut,
                scale_steps=scale_steps
            )

            C.analyze_metal_strain()
            C.analyze_porosity()
            # C.analyze_formation_energy()
            C.analyze_ligand_strain(
                # Assumes only one type of metal atom.
                metal_atom_no=[
                    cage_set.complex_dicts[i]['metal_atom_no']
                    for i in cage_set.complex_dicts
                ][0],
                expected_ligands=expected_ligands
            )
            cage_set.built_cage_properties[C.name] = {
                'pw_prop': C.pw_data,
                'op_prop': C.op_data,
                # 'form_energy': C.FE,
                'li_prop': C.ls_data,
                'fa_prop': C.fa_data,
            }
            # Dump to JSON.
            cage_set.dump_properties()

        cage_sets.append(cage_set)

    return cage_sets


def main():
    first_line = (
        'Usage: build_cage_library.py lig_lib_file prism_lib_file '
        'compl_lib_file lig_directory compl_directory'
    )
    if (not len(sys.argv) == 6):
        print(f"""
{first_line}

    lig_lib_file : (str)
        File containing ligand information (XXXXX)

    compl_lib_file : (str)
        File containing complex information (XXXXX).

    cage_lib_file : (str)
        File containing cage information (XXXXX).

    lig_directory : (str)
        Directory with required ligand structures.

    compl_directory : (str)
        Directory with required complex structures.

    """)
        sys.exit()
    else:
        lig_lib_file = sys.argv[1]
        compl_lib_file = sys.argv[2]
        cage_lib_file = sys.argv[3]
        ligand_directory = sys.argv[4]
        compl_directory = sys.argv[5]

    cage_lib = read_lib(cage_lib_file)
    compls = read_lib(compl_lib_file)
    ligs = read_lib(lig_lib_file)

    # Build and optimise all organic molecules in lib.
    cage_sets = build_cages(
        ligs,
        compls,
        cage_lib,
        ligand_directory,
        compl_directory
    )
    analyse_cages(cage_sets)


if __name__ == "__main__":
    main()
