#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to align crystal structures to calculated structures.

Author: Andrew Tarzia

Date Created: 03 Feb 2022

"""


import sys
import os
import numpy as np
from scipy.spatial.distance import cdist
from itertools import product

import stk
import spindry as spd

from utilities import read_lib
from xtalcage import XtalCage



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
            'cage_set': expt,
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
        'jd354': 1,
        'jd370': 1,
    }

    for xtal in xtals:
        print(f'---- doing: {xtal}')
        aid_ = aligned_ids[xtal]
        xtal_struct = stk.BuildingBlock.init_from_file(
            f'{xtal}_a_{aid_}_xtal.mol'
        )
        comp_struct = stk.BuildingBlock.init_from_file(
            f'{xtal}_a_{aid_}_comp.mol'
        )
        print(xtal_struct)
        print(comp_struct)
        raise SystemExit()


if __name__ == "__main__":
    main()
