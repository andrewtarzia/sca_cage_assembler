#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for new topology.

Author: Andrew Tarzia

Date Created: 28 May 2021

"""

import stk


class M8L6KnotPrism(stk.cage.Cage):
    """
    Represents a cage topology graph.

    Metal building blocks with three functional groups are
    required for this topology.

    Ligand building blocks with four functional groups are required for
    this topology.

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        # 3 Functional groups.
        stk.cage.vertices._NonLinearCageVertex(
            id=0,
            position=[-8.0, 5.5, 12.2],
            use_neighbor_placement=False,
        ),
        stk.cage.vertices._NonLinearCageVertex(
            id=1,
            position=[-0.2, 13.8, 2.9],
            use_neighbor_placement=False,
        ),
        stk.cage.vertices._NonLinearCageVertex(
            id=2,
            position=[13.1, 7.3, 4.1],
            use_neighbor_placement=False,
        ),
        stk.cage.vertices._NonLinearCageVertex(
            id=3,
            position=[-4.4, -10.3, 11.2],
            use_neighbor_placement=False,
        ),
        stk.cage.vertices._NonLinearCageVertex(
            id=4,  # 10,
            position=[-8.0, -5.5, -12.2],
            use_neighbor_placement=False,
        ),
        stk.cage.vertices._NonLinearCageVertex(
            id=5,  # 11,
            position=[-4.4, 10.3, -11.2],
            use_neighbor_placement=False,
        ),
        stk.cage.vertices._NonLinearCageVertex(
            id=6,  # 5,
            position=[13.1, -7.3, -4.1],
            use_neighbor_placement=False,
        ),
        stk.cage.vertices._NonLinearCageVertex(
            id=7,  # 6,
            position=[-0.2, -13.8, -2.9],
            use_neighbor_placement=False,
        ),

        # 4 Functional groups.
        stk.cage.vertices._NonLinearCageVertex(
            id=8,  # 4,
            position=[3.4, -1.4, 6.0],
            use_neighbor_placement=False,
        ),
        stk.cage.vertices._NonLinearCageVertex(
            id=9,  # 7,
            position=[2.1, -1.0, 3.7],
            use_neighbor_placement=False,
        ),
        stk.cage.vertices._NonLinearCageVertex(
            id=10,  # 8,
            position=[-7.1, 0.0, 0.0],
            use_neighbor_placement=False,
        ),
        stk.cage.vertices._NonLinearCageVertex(
            id=11,  # 9,
            position=[-4.4, -0.0, -0.0],
            use_neighbor_placement=False,
        ),
        stk.cage.vertices._NonLinearCageVertex(
            id=12,
            position=[3.4, 1.4, -6.0],
            use_neighbor_placement=False,
        ),
        stk.cage.vertices._NonLinearCageVertex(
            id=13,
            position=[2.1, 1.0, -3.7],
            use_neighbor_placement=False,
        ),
    )

    _edge_prototypes = (
        stk.Edge(
            id=0,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[8],
        ),
        stk.Edge(
            id=1,
            vertex1=_vertex_prototypes[3],
            vertex2=_vertex_prototypes[8],
        ),
        stk.Edge(
            id=2,
            vertex2=_vertex_prototypes[6],
            vertex1=_vertex_prototypes[8],
        ),
        stk.Edge(
            id=3,
            vertex1=_vertex_prototypes[2],
            vertex2=_vertex_prototypes[8],
        ),
        stk.Edge(
            id=4,
            vertex1=_vertex_prototypes[7],
            vertex2=_vertex_prototypes[9],
        ),
        stk.Edge(
            id=5,
            vertex1=_vertex_prototypes[3],
            vertex2=_vertex_prototypes[9],
        ),
        stk.Edge(
            id=6,
            vertex1=_vertex_prototypes[1],
            vertex2=_vertex_prototypes[9],
        ),
        stk.Edge(
            id=7,
            vertex1=_vertex_prototypes[2],
            vertex2=_vertex_prototypes[9],
        ),
        stk.Edge(
            id=8,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[10],
        ),
        stk.Edge(
            id=9,
            vertex1=_vertex_prototypes[3],
            vertex2=_vertex_prototypes[10],
        ),
        stk.Edge(
            id=10,
            vertex1=_vertex_prototypes[1],
            vertex2=_vertex_prototypes[11],
        ),
        stk.Edge(
            id=11,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[11],
        ),
        stk.Edge(
            id=12,
            vertex2=_vertex_prototypes[4],
            vertex1=_vertex_prototypes[10],
        ),
        stk.Edge(
            id=13,
            vertex2=_vertex_prototypes[5],
            vertex1=_vertex_prototypes[10],
        ),
        stk.Edge(
            id=14,
            vertex1=_vertex_prototypes[7],
            vertex2=_vertex_prototypes[11],
        ),
        stk.Edge(
            id=15,
            vertex2=_vertex_prototypes[4],
            vertex1=_vertex_prototypes[11],
        ),
        stk.Edge(
            id=16,
            vertex1=_vertex_prototypes[4],
            vertex2=_vertex_prototypes[12],
        ),
        stk.Edge(
            id=17,
            vertex1=_vertex_prototypes[5],
            vertex2=_vertex_prototypes[12],
        ),
        stk.Edge(
            id=18,
            vertex1=_vertex_prototypes[2],
            vertex2=_vertex_prototypes[12],
        ),
        stk.Edge(
            id=19,
            vertex1=_vertex_prototypes[6],
            vertex2=_vertex_prototypes[12],
        ),
        stk.Edge(
            id=20,
            vertex1=_vertex_prototypes[1],
            vertex2=_vertex_prototypes[13],
        ),
        stk.Edge(
            id=21,
            vertex1=_vertex_prototypes[5],
            vertex2=_vertex_prototypes[13],
        ),
        stk.Edge(
            id=22,
            vertex1=_vertex_prototypes[7],
            vertex2=_vertex_prototypes[13],
        ),
        stk.Edge(
            id=23,
            vertex1=_vertex_prototypes[6],
            vertex2=_vertex_prototypes[13],
        ),
    )

    def _get_scale(self, building_block_vertices):
        return 2
