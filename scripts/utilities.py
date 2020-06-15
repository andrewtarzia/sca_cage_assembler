#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Utilities module.

Author: Andrew Tarzia

Date Created: 15 Mar 2020
"""

import json
from itertools import combinations

from atools import get_atom_distance


def read_lib(lib_file):
    """
    Read lib file.

    Returns dictionary.

    """

    print(f'reading {lib_file}')
    with open(lib_file, 'rb') as f:
        lib = json.load(f)

    return lib


def calculate_binding_AR(mol):
    """
    Calculate ligand aspect ratio based on binder positions.

    Defined as:

    """

    if mol.get_num_functional_groups() != 4:
        return None

    binder_atom_ids = [
        list(fg.get_bonder_ids())
        for fg in mol.get_functional_groups()
    ]
    binder_atom_dists = sorted(
        [
            (idx1, idx2, get_atom_distance(mol, idx1, idx2))
            for idx1, idx2 in combinations(binder_atom_ids, r=2)
        ],
        key=lambda a: a[2]
    )

    far_binder_pair = (
        binder_atom_dists[-1][0],
        binder_atom_dists[-1][1]
    )
    ARs = []
    for fg_id in far_binder_pair:
        ds = sorted([
            i[2] for i in binder_atom_dists
            if fg_id in (i[0], i[1])
        ])
        AR = ds[1]/min(ds)
        ARs.append(AR)

    ligand_AR = sum(ARs)/len(ARs)
    return ligand_AR
