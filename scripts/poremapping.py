#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run poremapper on aligned structures.

Author: Andrew Tarzia

Date Created: 17 Feb 2022

"""

import sys

import stk
import pore_mapper as pm

from utilities import read_lib


def run_poremapper(structure, file_prefix):
    structure.write(f'{file_prefix}_host.xyz')
    # Read in host from xyz file.
    host = pm.Host.init_from_xyz_file(path=f'{file_prefix}_host.xyz')
    host = host.with_centroid([0., 0., 0.])

    # Define calculator object.
    calculator = pm.Inflater(bead_sigma=1.2)
    final_result = calculator.get_inflated_blob(host=host)
    pore = final_result.pore

    # Do final structure.
    host.write_xyz_file(f'{file_prefix}_host.xyz')
    pore.write_xyz_file(f'{file_prefix}_pore.xyz')
    return pore.get_volume()


def main():
    first_line = (
        'Usage: analyse_crystal_structures.py '
        'expt_lib_file'
    )
    if (not len(sys.argv) == 2):
        print(f"""
{first_line}

    expt_lib_file : (str)
        File containing experimental symmetry  information (XXXXX).

    """)
        sys.exit()
    else:
        expt_lib_file = sys.argv[1]

    expt_data = read_lib(expt_lib_file)

    # List of the xtal structures and their corresponding names.
    xtals = {}
    for expt in expt_data:
        xtals[expt_data[expt]['xtal_struct_name']] = {
            'cage_set': expt_data[expt]['cage_set'],
            'symmetry_name': expt_data[expt]['symmetry'],
            'ligand_name': expt_data[expt]['ligand_name'],
            'complexes': tuple(expt_data[expt]['complexes']),
        }

    # Hardcode the best aligned ID.
    aligned_ids = {
        'jd235': 0,
        'jd257': 14,
        'jd301': 6,
        'jd326': 0,
        'jd354': 15,
        'jd370l': 0,
        'jd370d': 9,
    }

    for xtal in xtals:
        aid_ = aligned_ids[xtal]
        xtal_struct = stk.BuildingBlock.init_from_file(
            f'{xtal}_a_{aid_}_xtal.mol'
        )
        xtal_porevolume = run_poremapper(
            structure=xtal_struct,
            file_prefix=f'{xtal}_pm_xtal',
        )

        comp_struct = stk.BuildingBlock.init_from_file(
            f'{xtal}_a_{aid_}_comp.mol'
        )
        comp_porevolume = run_poremapper(
            structure=comp_struct,
            file_prefix=f'{xtal}_pm_comp',
        )

        print(f'---- {xtal}: {xtal_porevolume} vs. {comp_porevolume}')


if __name__ == "__main__":
    main()
