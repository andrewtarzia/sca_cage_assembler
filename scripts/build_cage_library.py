#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build complex library.

Author: Andrew Tarzia

Date Created: 27 Jan 2020
"""

import sys
from itertools import combinations
from os.path import exists
import numpy as np
import matplotlib.pyplot as plt

from cage_set import HoCube, HetPrism
from utilities import read_lib, heatmap, annotate_heatmap


def plot_set_heatmap_vs_aniso(data, Y_name, ylabel, ylim, filename):

    symmetries = [
        'o1', 'th1', 'th2', 't1', 's61', 's62', 'd31', 'd32',
        'c2v', 'c2h',
    ]
    col_names = {i: j['aspect_ratio'] for i, j in data.items()}
    col_names = {
        i: j for i, j in sorted(
            col_names.items(), key=lambda item: item[1]
        )
    }
    aspect_ratios = sorted(col_names.values())
    out_values = np.zeros((len(symmetries), len(aspect_ratios)))

    fig, ax = plt.subplots()

    for name in data:
        col = list(col_names.keys()).index(name)
        # Iterate over all cages in set.
        for cage in data[name][Y_name]:
            Y = data[name][Y_name][cage]
            symm = cage.split('_')[-1]
            row = symmetries.index(symm)

            # Assign to zeros matrix.
            if out_values[row][col] != 0:
                raise ValueError()

            out_values[row][col] = Y

    im, cbar = heatmap(
        data=out_values,
        row_labels=symmetries,
        col_labels=[f'{round(i, 2)}' for i in col_names.values()],
        ax=ax,
        cmap='Blues',
        cbarlabel=ylabel
    )
    annotate_heatmap(im, valfmt="{x:.1f}", fontsize=8)
    ax.set_xlabel('ligand aspect ratio [1:X]')
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
        'lse_sum': {},
        'min_imine_torsions': {},
        'max_ligand_distortion': {},
        'max_diff_face_aniso': {},
        'max_ML_length': {},
        'formation_energies': {},
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
        measures['lse_sum'][C.name] = sum([
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
            100*((i[2] - i[3]) / i[2])
            for i in C_data['fa_prop']
        ])
        measures['max_ML_length'][C.name] = max([
            i for i in C_data['bl_prop']
        ])
        measures['formation_energies'][C.name] = C_data['fe_prop']

    print('minimum OPs:', measures['oct_op'])
    print('sum LSE:', measures['lse_sum'])
    print('min imine torsion:', measures['min_imine_torsions'])
    print('max core planarities:', measures['max_ligand_distortion'])
    print('max diff in face aniso:', measures['max_diff_face_aniso'])
    print('max metal-ligand distance:', measures['max_ML_length'])
    print('formation energies:', measures['formation_energies'])
    print('----------------------------------------------------')

    plottables = {
        'oct_op': {
            'data': measures['oct_op'],
            'ylabel': r'min. $q_{\mathrm{oct}}$',
            'ylim': (0, 1),
            'filename': f'{cage_set.name}_minOPs.pdf'
        },
        'lse_sum': {
            'data': {
                i: (
                    measures['lse_sum'][i]
                    - min(measures['lse_sum'].values())
                )
                for i in measures['lse_sum']
            },
            'ylabel': r'rel. sum strain energy [kJ/mol]',
            'ylim': (-4, 500),
            'filename': f'{cage_set.name}_sumLSE.pdf'
        },
        'min_imine_torsions': {
            'data': measures['min_imine_torsions'],
            'ylabel': r'min. imine torsion [degrees]',
            'ylim': (0, 185),
            'filename': f'{cage_set.name}_mintors.pdf'
        },
        'max_ligand_distortion': {
            'data': measures['max_ligand_distortion'],
            'ylabel': r'max. ligand distortion [$\mathrm{\AA}$]',
            'ylim': (0, 185),
            'filename': f'{cage_set.name}_maxdistortion.pdf'
        },
        'max_diff_face_aniso': {
            'data': measures['max_diff_face_aniso'],
            'ylabel': r'max. $\Delta$opposing face anisotropy [%]',
            'ylim': (-10, 100),
            'filename': f'{cage_set.name}_maxfadiff.pdf'
        },
        'max_ML_length': {
            'data': measures['max_ML_length'],
            'ylabel': r'max. N-Zn bond length [$\mathrm{\AA}$]',
            'ylim': (2, 3),
            'filename': f'{cage_set.name}_maxmld.pdf'
        },
        'relfe': {
            'data': {
                i: (
                    measures['formation_energies'][i]
                    - min(measures['formation_energies'].values())
                )
                for i in measures['formation_energies']
            },
            'ylabel': r'rel. formation energy [kJ/mol]',
            'ylim': (-10, 1000),
            'filename': f'{cage_set.name}_relfe.pdf'
        },
    }

    for p1, p2 in combinations(plottables, 2):
        p1_dict = plottables[p1]
        if not exists(p1_dict['filename']):
            cage_set.plot_Y(
                data=p1_dict['data'],
                ylabel=p1_dict['ylabel'],
                ylim=p1_dict['ylim'],
                filename=p1_dict['filename']
            )
        p2_dict = plottables[p2]
        if not exists(p2_dict['filename']):
            cage_set.plot_Y(
                data=p2_dict['data'],
                ylabel=p2_dict['ylabel'],
                ylim=p2_dict['ylim'],
                filename=p2_dict['filename']
            )

        p1_p2_filename = f'{cage_set.name}_{p1}_{p2}.pdf'
        if not exists(p1_p2_filename):
            cage_set.plot_Y_C(
                data=p1_dict['data'],
                ylabel=p1_dict['ylabel'],
                ylim=p1_dict['ylim'],
                data_C=p2_dict['data'],
                clabel=p2_dict['ylabel'],
                clim=p2_dict['ylim'],
                filename=p1_p2_filename
            )
        p2_p1_filename = f'{cage_set.name}_{p2}_{p1}.pdf'
        if not exists(p2_p1_filename):
            cage_set.plot_Y_C(
                data=p2_dict['data'],
                ylabel=p2_dict['ylabel'],
                ylim=p2_dict['ylim'],
                data_C=p1_dict['data'],
                clabel=p1_dict['ylabel'],
                clim=p1_dict['ylim'],
                filename=p2_p1_filename
            )

    return measures


def analyse_cages(cage_sets):

    AR_data = {}
    for cage_set in cage_sets:
        if isinstance(cage_set, HetPrism):
            het_prism_analysis(cage_set)
        if isinstance(cage_set, HoCube):
            measures = homo_cube_analysis(cage_set)
            # For all HoCube sets, collect ligand aspect ratio data.
            AR_data[cage_set.name] = {
                'min_OPs': measures['oct_op'],
                'lse_sum': {
                    i: (
                        measures['lse_sum'][i] -
                        min(measures['lse_sum'].values())
                    )
                    for i in measures['lse_sum']
                },
                'min_tor': measures['min_imine_torsions'],
                'max_dis': measures['max_ligand_distortion'],
                'aspect_ratio': cage_set.ligand_aspect_ratio,
                'max_rfa': measures['max_diff_face_aniso'],
                'max_mld': measures['max_ML_length'],
                'rel_fe': {
                    i: (
                        measures['formation_energies'][i]
                        - min(measures['formation_energies'].values())
                    )
                    for i in measures['formation_energies']
                },
            }

    # Plot ligand aspect ratio data.
    tests = {
        # Test: (ylabel, ylim)
        'min_OPs': (r'min. $q_{\mathrm{oct}}$', (0, 1)),
        'lse_sum': (r'rel. sum strain energy [kJ/mol]', (-4, 500)),
        'min_tor': (r'min. imine torsion [degrees]', (0, 185)),
        'max_dis': (
            r'max. ligand distortion [$\mathrm{\AA}$]', (0, 185)
        ),
        'max_rfa': (
            r'max. $\Delta$opposing face anisotropy [%]', (-10, 100)
        ),
        'max_mld': (
            r'max. N-Zn bond length [$\mathrm{\AA}$]', (2, 3)
        ),
        'rel_fe': (r'rel. formation energy [kJ/mol]', (-10, 1000)),
    }
    if len(AR_data) > 0:
        for t in tests:
            plot_set_heatmap_vs_aniso(
                data=AR_data,
                Y_name=t,
                ylabel=tests[t][0],
                ylim=tests[t][1],
                filename=f'plotset_{t}_VA.pdf'
            )


def build_cages(
    ligands,
    complexes,
    cage_set_lib,
    ligand_directory,
    complex_directory,
    read_data,
):

    cage_sets = []
    for name in cage_set_lib:
        cage_set_c = cage_set_lib[name]
        compl_names = cage_set_c['corners']
        comps = {i: complexes[i] for i in compl_names}

        if cage_set_c['heteroleptic']:
            cage_set = HetPrism(
                name=name,
                cage_set_dict=cage_set_c,
                complex_dicts=comps,
                ligand_dicts=ligands,
                ligand_dir=ligand_directory,
                complex_dir=complex_directory
            )
        else:
            cage_set = HoCube(
                name=name,
                cage_set_dict=cage_set_c,
                complex_dicts=comps,
                ligand_dicts=ligands,
                ligand_dir=ligand_directory,
                complex_dir=complex_directory
            )

        if read_data and exists(cage_set.properties_file):
            cage_set.load_properties()
        else:
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
                C.analyze_ligand_strain(
                    # Assumes only one type of metal atom.
                    metal_atom_no=[
                        cage_set.complex_dicts[i]['metal_atom_no']
                        for i in cage_set.complex_dicts
                    ][0],
                    expected_ligands=expected_ligands,
                    free_e=default_free_e,
                )
                cage_set.built_cage_properties[C.name] = {
                    'pw_prop': C.pw_data,
                    'op_prop': C.op_data,
                    'fe_prop': C.fe_data,
                    'li_prop': C.ls_data,
                    'fa_prop': C.fa_data,
                    'bl_prop': C.bl_data
                }
                # Dump to JSON.
                cage_set.dump_properties()

        cage_sets.append(cage_set)

    return cage_sets


def main():
    first_line = (
        'Usage: build_cage_library.py lig_lib_file prism_lib_file '
        'compl_lib_file lig_directory compl_directory read_data'
    )
    if (not len(sys.argv) == 7):
        print(f"""
{first_line}

    lig_lib_file : (str)
        File containing ligand information (XXXXX)

    compl_lib_file : (str)
        File containing complex information (XXXXX).

    cage_set_lib_file : (str)
        File containing cage information (XXXXX).

    lig_directory : (str)
        Directory with required ligand structures.

    compl_directory : (str)
        Directory with required complex structures.

    read_data : (str)
        't' if cage analysis can be read from CageSet.properties_file.
        All other strings gives False.

    """)
        sys.exit()
    else:
        lig_lib_file = sys.argv[1]
        compl_lib_file = sys.argv[2]
        cage_set_lib_file = sys.argv[3]
        ligand_directory = sys.argv[4]
        compl_directory = sys.argv[5]
        read_data = True if sys.argv[6] == 't' else False

    cage_set_lib = read_lib(cage_set_lib_file)
    compls = read_lib(compl_lib_file)
    ligs = read_lib(lig_lib_file)

    # Build and optimise all organic molecules in lib.
    cage_sets = build_cages(
        ligs,
        compls,
        cage_set_lib,
        ligand_directory,
        compl_directory,
        read_data,
    )
    analyse_cages(cage_sets)


if __name__ == "__main__":
    main()
