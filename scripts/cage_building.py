#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Modules defining and building the Cage class and subclasses.

Author: Andrew Tarzia

Date Created: 03 Mar 2020

"""

import stk

from knot_topology import M8L6KnotPrism


def m4l4spacer_graph(building_blocks, vertex_alignments):

    topology_graph = stk.cage.M4L4Tetrahedron(
        building_blocks=building_blocks,
        vertex_alignments=vertex_alignments,
        num_processes=2,
    )

    return topology_graph


def m8l6_graph(building_blocks, vertex_alignments):

    topology_graph = stk.cage.M8L6Cube(
        building_blocks=building_blocks,
        vertex_alignments=vertex_alignments,
        num_processes=2,
    )

    return topology_graph


def m8l6knot_graph(building_blocks, vertex_alignments):

    topology_graph = M8L6KnotPrism(
        building_blocks=building_blocks,
        vertex_alignments=vertex_alignments,
        num_processes=2,
        optimizer=stk.MCHammer(),
    )

    return topology_graph


def m6l2l3_graph(building_blocks, vertex_alignments):

    topology_graph = stk.cage.M6L2L3Prism(
        building_blocks=building_blocks,
        vertex_alignments=vertex_alignments,
        num_processes=2,
    )

    return topology_graph


def available_topologies(string):
    """
    Get stk function of desired topology.

    """

    topologies = {
        'm4l4spacer': m4l4spacer_graph,
        'm8l6face': m8l6_graph,
        'm6l2l3': m6l2l3_graph,
        'm8l6knot': m8l6knot_graph,
    }

    try:
        return topologies[string]
    except KeyError:
        raise KeyError(f'{string} not in {topologies.keys()}')
