#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Modules defining and building the Cage class and subclasses.

Author: Andrew Tarzia

Date Created: 03 Mar 2020

"""

import stk


def m4l4spacer_graph(metal_corner, ligand, ligand2=None):

    raise NotImplementedError('see m8l6 for example')

    if ligand.get_num_functional_groups() != 3:
        raise ValueError(f'{ligand} does not have 3 functional groups')

    cage = stk.ConstructedMolecule(
        topology_graph=stk.cage.M4L4Tetrahedron(
            building_blocks={
                metal_corner: range(4),
                ligand: range(4, 8)
            },
        )
    )

    return cage


def m8l6_graph(building_blocks, vertex_alignments):

    topology_graph = stk.cage.M8L6Cube(
        building_blocks=building_blocks,
        vertex_alignments=vertex_alignments,
        num_processes=2,
    )

    return topology_graph


def m6l2l3_graph(metal_corner, ligand, ligand2):

    raise NotImplementedError('see m8l6 for example')

    lig_num_fgs = ligand.get_num_functional_groups()
    lig2_num_fgs = ligand2.get_num_functional_groups()

    if lig_num_fgs == 3 and lig2_num_fgs == 4:
        tri_ligand = ligand
        tet_ligand = ligand2
    elif lig_num_fgs == 4 and lig2_num_fgs == 3:
        tri_ligand = ligand2
        tet_ligand = ligand
    else:
        raise ValueError(
            f'{ligand} or {ligand2} does not have 3 or 4 functional '
            'groups'
        )

    cage = stk.ConstructedMolecule(
        topology_graph=stk.cage.M6L2L3Prism(
            building_blocks={
                metal_corner: range(6),
                tri_ligand: range(6, 8),
                tet_ligand: range(8, 11)
            },
        )
    )

    return cage


def available_topologies(string):
    """
    Get stk function of desired topology.

    """

    topologies = {
        'm4l4spacer': m4l4spacer_graph,
        'm8l6face': m8l6_graph,
        'm6l2l3': m6l2l3_graph,
    }

    try:
        return topologies[string]
    except KeyError:
        raise KeyError(f'{string} not in {topologies.keys()}')


def defined_face_sets(string):
    """
    Define which vertices in a topology share a face.

    """

    topologies = {
        'm8l6face': {
            '001': {
                'vertices': (0, 1, 2, 3),
                'connected': (
                    (0, 1), (1, 2), (2, 3), (3, 0)
                ),
                'opposite': '00-1'
            },
            '100': {
                'vertices': (0, 1, 4, 5),
                'connected': (
                    (0, 1), (1, 5), (5, 4), (4, 0)
                ),
                'opposite': '-100'
            },
            '010': {
                'vertices': (0, 3, 4, 7),
                'connected': (
                    (0, 3), (3, 7), (7, 4), (4, 0)
                ),
                'opposite': '0-10'
            },
            '00-1': {
                'vertices': (4, 5, 6, 7),
                'connected': (
                    (4, 5), (5, 6), (6, 7), (7, 4)
                ),
                'opposite': '001'
            },
            '-100': {
                'vertices': (2, 3, 6, 7),
                'connected': (
                    (2, 3), (3, 7), (7, 6), (6, 2)
                ),
                'opposite': '100'
            },
            '0-10': {
                'vertices': (1, 2, 5, 6),
                'connected': (
                    (1, 2), (2, 6), (6, 5), (5, 1)
                ),
                'opposite': '010'
            },
        },
    }

    try:
        return topologies[string]
    except KeyError:
        raise KeyError(f'{string} not in {topologies.keys()}')
