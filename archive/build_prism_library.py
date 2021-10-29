#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build HetPrism library.

Author: Andrew Tarzia

Date Created: 27 Jan 2020
"""

import sys
from os.path import exists

from cage_set import HetPrism
from cage_analysis import analyse_cages
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
        cage_set_c = cage_set_lib[name]
        compl_names = cage_set_c['corners']
        comps = {i: complexes[i] for i in compl_names}

        cage_set = HetPrism(
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
                    distance_cut = 3.0
                    scale_steps = True
                    expected_ligands = 2
                elif C.topology_string == 'm8l6face':
                    distance_cut = 2.5
                    scale_steps = False
                    expected_ligands = 1
                else:
                    distance_cut = 2.0
                    scale_steps = True
                    expected_ligands = 1
                target_bond_length = 1.2
                num_steps = 2000
                step_size = 0.25

                C.optimize(
                    free_e=default_free_e,
                    step_size=step_size,
                    target_bond_length=target_bond_length,
                    num_steps=num_steps
                )

                continue

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
        'Usage: build_prism_library.py lig_lib_file prism_lib_file '
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

    # Build and optimise all organic molecules in lib.
    cage_sets = build_cages(
        ligands=ligands,
        complexes=complexes,
        cage_set_lib=cage_set_lib,
        ligand_directory=ligand_directory,
        complex_directory=complex_directory,
        read_data=read_data,
    )
    analyse_cages(cage_sets)


if __name__ == "__main__":
    main()