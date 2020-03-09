#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build complex library.

Author: Andrew Tarzia

Date Created: 27 Jan 2020
"""

import sys
import json
import numpy as np

import cage_building


def read_lib(lib_file):
    """
    Read lib file.

    Returns dictionary.

    """

    with open(lib_file, 'r') as f:
        lib = json.load(f)

    return lib


def analyse_cages(het_cages):

    for het_cage in het_cages:
        het_cage.load_properties()
        # Compare average pore volume of each of the three topologies.
        three_top = {'m4l4spacer': [], 'm8l6face': [], 'm6l2l3': []}
        built_prop = het_cage.built_cage_properties
        all_pore_volumes = []
        all_min_OPs = []
        all_topo_strs = []
        for C in het_cage.cages_to_build:
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
        for C in het_cage.cages_to_build:
            C_data = built_prop[C.name]
            TOPO = C.topology_string
            # Get minimium octahedral OP of the metal that is in the
            # complex building block.
            print(het_cage.complex_dicts)
            atom_no_of_interest = list(set([
                int(het_cage.complex_dicts[i]['metal_atom_no'])
                for i in het_cage.complex_dicts
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
        het_cage.plot_min_OPs_avg_PV(
            X=all_pore_volumes,
            Y=all_min_OPs,
            T=all_topo_strs
        )
        sys.exit()


def build_cages(
    complexes,
    prisms,
    ligand_directory,
    complex_directory
):

    cages = []
    for name in prisms:
        pris = prisms[name]
        compl_names = pris['corners'].strip(')(').split(', ')
        print(compl_names)
        comps = {i: complexes[i] for i in compl_names}
        print(pris)
        print(comps)

        het_cage = cage_building.HetPrism(
            name=name,
            prism_dict=pris,
            complex_dicts=comps,
            ligand_dir=ligand_directory,
            complex_dir=complex_directory
        )
        print('tob,..............', het_cage.cages_to_build)
        for C in het_cage.cages_to_build:
            print(C)
            C.build()
            default_free_e = C.free_electron_options[0]
            print(C.free_electron_options, default_free_e)
            if 'm6l2l3' in C.name:
                # Use a slightly different collapser threshold for
                # prism.
                step_size = 0.05
                distance_cut = 3.0
            else:
                step_size = 0.05
                distance_cut = 2.0
            C.optimize(
                free_e=default_free_e,
                step_size=step_size,
                distance_cut=distance_cut
            )
            C.analyze_cage_geometry()
            C.analyze_cage_porosity()
            het_cage.built_cage_properties[C.name] = {
                'pw_prop': C.pw_data,
                'op_prop': C.op_data,
            }
            # Dump to JSON.
            het_cage.dump_properties()

        cages.append(het_cage)

    return cages


def main():
    first_line = (
        'Usage: build_cage_library.py prism_lib_file compl_lib_file'
        'lig_directory compl_directory'
    )
    if (not len(sys.argv) == 5):
        print(f"""
{first_line}

    compl_lib_file : (str)
        File containing complex information (XXXXX).

    prism_lib_file : (str)
        File containing prism information (XXXXX).

    lig_directory : (str)
        Directory with required ligand structures.

    compl_directory : (str)
        Directory with required complex structures.

    """)
        sys.exit()
    else:
        compl_lib_file = sys.argv[1]
        prism_lib_file = sys.argv[2]
        ligand_directory = sys.argv[3]
        compl_directory = sys.argv[4]

    print(f'reading {prism_lib_file}')
    prisms = read_lib(prism_lib_file)
    print(f'reading {compl_lib_file}')
    compls = read_lib(compl_lib_file)

    # Build and optimise all organic molecules in lib.
    cages = build_cages(
        compls,
        prisms,
        ligand_directory,
        compl_directory
    )

    analyse_cages(cages)

    sys.exit()


if __name__ == "__main__":
    main()
