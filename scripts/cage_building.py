#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Modules defining and building the Cage class and subclasses.

Author: Andrew Tarzia

Date Created: 03 Mar 2020

"""

import stk


def m8l6_graph(building_blocks, vertex_alignments):

    topology_graph = stk.cage.M8L6Cube(
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
        'm8l6face': m8l6_graph,
    }

    try:
        return topologies[string]
    except KeyError:
        raise KeyError(f'{string} not in {topologies.keys()}')
