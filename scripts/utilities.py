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
        Average ratio of the two shortest binding atom-binding atom
        distances eminating from each binding atom in the molecule.

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


def calculate_paired_face_anisotropies(mol, metal_atom_ids, face_sets):
    """
    Calculate face anisotropy of opposing sides of a prism.

    Defined as:
        Ratio of metal-metal distances defined by two edges eminating
        from a chosen metal on each face.

    """

    # Get all metal-metal distances with metal atom ids.
    metal_atom_dists = sorted(
        [
            (idx1, idx2, get_atom_distance(mol, idx1, idx2))
            for idx1, idx2 in combinations(metal_atom_ids, r=2)
        ],
        key=lambda a: a[2]
    )

    face_anisotropies = {}
    for fs in face_sets:
        fsv = face_sets[fs]['vertices']
        fsc = face_sets[fs]['connected']
        fs_atom_ids = tuple(metal_atom_ids[i] for i in fsv)
        fs_distances = tuple(
            i for i in metal_atom_dists
            if i[0] in fs_atom_ids and i[1] in fs_atom_ids
        )
        # Pick one metal.
        idx = fsv[0]
        conn = [
            j
            for i in fsc if idx in i
            for j in i if j != idx
        ]
        fs_idx = metal_atom_ids[idx]
        fs_conn = [metal_atom_ids[i] for i in conn]

        # Define aniso based on the distances between its two
        # connections.
        pair1 = (fs_idx, fs_conn[0])
        pair2 = (fs_idx, fs_conn[1])
        d1, = (
            i[2]
            for i in fs_distances
            if i[0] in pair1 and i[1] in pair1
        )
        d2, = (
            i[2]
            for i in fs_distances
            if i[0] in pair2 and i[1] in pair2
        )
        face_aniso = d2 / d1 if d2 > d1 else d1 / d2
        face_anisotropies[fs] = face_aniso

    # Pair the face anisotropies of opposing faces.
    paired_face_anisotropies = [
        (i, j, face_anisotropies[i], face_anisotropies[j])
        for i, j in combinations(face_anisotropies, r=2)
        if face_sets[i]['opposite'] == j
    ]

    return paired_face_anisotropies
