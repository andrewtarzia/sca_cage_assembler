#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build HoCube library.

Author: Andrew Tarzia

Date Created: 27 Jan 2020

"""

import sys
import logging
import os

from cage_analysis import write_csv

from env_set import read_envset_json, read_json

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)


def build_cages(
    ligands,
    complexes,
    cage_set_lib,
    environment,
):

    cages = []
    for name in cage_set_lib:
        logging.info(f'doing cage set: {name}')
        cage_set_definition = cage_set_lib[name]
        complexes = {i: complexes[i] for i in cage_set_definition['corners']}
        print(cage_set_definition)
        print(complexes)
        raise SystemExit()
        cage_set = HoCube(
            name=name,
            cage_set_dict=cage_set_c,
            complex_dicts=comps,
            ligand_dicts=ligands,
            ligand_dir=ligand_directory,
            complex_dir=complex_directory
        )

        if exists(cage_set.properties_file):
            cage_set.load_properties()

        if not read_data:
            for C in cage_set.cages_to_build:
                C.build()

                default_free_e = C.free_electron_options[0]
                # Use a slightly different collapser threshold for
                # different topologies.
                if C.topology_string in ['m8l6face', 'm8l6knot']:
                    step_size = 0.05
                    distance_cut = 2.5
                    scale_steps = False
                    expected_ligands = 1
                    # target_bond_length = 1.2
                    # num_steps = 2000
                    # step_size = 0.25
                else:
                    step_size = 0.05
                    distance_cut = 2.5
                    scale_steps = False
                    expected_ligands = 1
                C.save_bb_vector_xyzs(f'{C.unopt_file}.mol')

                # Check if structure has previously had optimisation
                # attempted with failure.
                try:
                    if (
                        cage_set.built_cage_properties[C.name][
                            'optimized'
                        ]
                        is False
                    ):
                        C.optimized = False
                        print(
                            'Warning!: xTB optimisation of '
                            f'{C.name} did not converge and structure '
                            'is being kept unoptimized.'
                        )
                except KeyError:
                    pass

                if C.optimized is None:
                    # Run optimisation - which handles successfully
                    # completed molecules.
                    C.optimize(
                        free_e=default_free_e,
                        step_size=step_size,
                        distance_cut=distance_cut,
                        scale_steps=scale_steps,
                    )
                if C.optimized is False:
                    cage_set.built_cage_properties[C.name] = {
                        'optimized': C.optimized,
                    }
                else:
                    C.save_bb_vector_xyzs(f'{C.opt_file}.mol')
                    C.analyze_ligand_strain(
                        # Assumes only one type of metal atom.
                        metal_atom_no=[
                            cage_set.complex_dicts[i]['metal_atom_no']
                            for i in cage_set.complex_dicts
                        ][0],
                        expected_ligands=expected_ligands,
                        free_e=default_free_e,
                    )

                    C.analyze_cube_likeness()
                    C.analyze_metal_strain()
                    C.analyze_porosity()
                    cage_set.built_cage_properties[C.name] = {
                        'optimized': C.optimized,
                        'pw_prop': C.pw_data,
                        'op_prop': C.op_data,
                        'fe_prop': C.fe_data,
                        'li_prop': C.ls_data,
                        'bl_prop': C.bl_data,
                        'cl_prop': C.cl_data,
                        'c8_prop': C.analyze_cube_shape(),
                    }
                # Dump to JSON.
                cage_set.dump_properties()

        cage_sets.append(cage_set)

    return cages


def main():
    if (not len(sys.argv) == 2):
        print("""
Usage: build_cage_library.py env_set_file

    env_set_file : (str)
        File containing environment set file.


    """)
        sys.exit()
    else:
        env_set_file = sys.argv[1]

    environment = read_envset_json(env_set_file)

    cage_set_lib = read_json(environment["cage_file"])
    complexes = read_json(environment["complex_file"])
    ligands = read_json(environment["ligand_file"])
    expt_data = read_json(environment['experiment_file'])

    # List of the cages that are known to form, including symmetry.
    experimentals = [
        f"C_{expt_data[i]['cage_set']}_{expt_data[i]['symmetry']}"
        for i in expt_data
    ]

    # Build and optimise all organic molecules in lib.
    cages = build_cages(
        ligands=ligands,
        complexes=complexes,
        cage_set_lib=cage_set_lib,
        environment=environment,
    )
    write_csv(cages, experimentals)


if __name__ == "__main__":
    main()
