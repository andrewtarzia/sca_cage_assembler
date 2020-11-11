#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build HoCube library.

Author: Andrew Tarzia

Date Created: 27 Jan 2020

"""

import sys
from os.path import exists

from cage_set import HoCube
from cage_analysis import analyse_cages, analyse_cage_sets
from utilities import read_lib


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
        print(f'------------- doing {name}:')
        cage_set_c = cage_set_lib[name]
        compl_names = cage_set_c['corners']
        comps = {i: complexes[i] for i in compl_names}

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
                if C.topology_string == 'm8l6face':
                    step_size = 0.05
                    distance_cut = 2.5
                    scale_steps = False
                    expected_ligands = 1
                    # target_bond_length = 1.2
                    # num_steps = 2000
                    # step_size = 0.25
                C.optimize(
                    free_e=default_free_e,
                    step_size=step_size,
                    distance_cut=distance_cut,
                    scale_steps=scale_steps,
                )
                C.analyze_cube_likeness()
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
                    'bl_prop': C.bl_data,
                    'cl_prop': C.cl_data,
                }
                # Dump to JSON.
                cage_set.dump_properties()

        cage_sets.append(cage_set)

    return cage_sets


def main():
    first_line = (
        'Usage: build_cube_library.py lig_lib_file prism_lib_file '
        'compl_lib_file lig_directory compl_directory read_data'
    )
    if (not len(sys.argv) == 7):
        print(f"""
{first_line}

    ligand_lib_file : (str)
        File containing ligand information (XXXXX)

    complex_lib_file : (str)
        File containing complex information (XXXXX).

    cage_set_lib_file : (str)
        File containing cage information (XXXXX).

    ligand_directory : (str)
        Directory with required ligand structures.

    complex_directory : (str)
        Directory with required complex structures.

    read_data : (str)
        't' if cage analysis can be read from CageSet.properties_file.
        All other strings gives False.

    """)
        sys.exit()
    else:
        ligand_lib_file = sys.argv[1]
        complex_lib_file = sys.argv[2]
        cage_set_lib_file = sys.argv[3]
        ligand_directory = sys.argv[4]
        complex_directory = sys.argv[5]
        read_data = True if sys.argv[6] == 't' else False

    cage_set_lib = read_lib(cage_set_lib_file)
    complexes = read_lib(complex_lib_file)
    ligands = read_lib(ligand_lib_file)

    # List of the cages that are known to form, including symmetry.
    experimentals = [
        'C_cl1_quad2_12_th1',
        'C_cl1_quad2_12_th2',
        'C_cl1_quad2_16_th1',
        'C_cl1_quad2_16_th2',
        'C_cl1_quad2_3_th1',
        'C_cl1_quad2_3_th2',
        'C_cl1_quad2_8_s61',
        'C_cl1_quad2_8_s62',
    ]

    # Build and optimise all organic molecules in lib.
    cage_sets = build_cages(
        ligands=ligands,
        complexes=complexes,
        cage_set_lib=cage_set_lib,
        ligand_directory=ligand_directory,
        complex_directory=complex_directory,
        read_data=read_data,
    )
    analyse_cage_sets(cage_sets)
    analyse_cages(cage_sets, experimentals)


if __name__ == "__main__":
    main()
