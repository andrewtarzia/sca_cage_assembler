#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Modules defining and building the Cage class and subclasses.

Author: Andrew Tarzia

Date Created: 03 Mar 2020

"""

from itertools import product
from os.path import exists, join

import stk


def available_topologies(string):
    """
    Get stk function of desired topology.

    """

    topologies = {
        'm4l4spacer': stk.cage.M4L4_Oct_Spacer(),
        'm8l6face': stk.cage.M8L6_Oct_Face(),
        'm6l2l3': stk.cage.M6L2L3_Oct()
    }

    try:
        return string, topologies[string]
    except KeyError:
        raise KeyError(f'{string} not in {topologies.keys()}')


class Cage:
    """
    Generic class that builds and analyses stk.ConstructuedMolecules.

    """

    def __init__(self, name, bbs, topology, bb_vertices):
        self.name = name
        self.bbs = bbs
        self.topology = topology
        self.bb_vertices = bb_vertices

    def build(self):
        print(f'....building {self.name}')
        cage = stk.ConstructedMolecule(
            building_blocks=self.bbs,
            topology_graph=self.topology,
            building_block_vertices=self.bb_vertices
        )
        cage.write(f'{self.name}_unopt.mol')
        cage.dump(f'{self.name}_unopt.json')
        self.cage = cage

    def optimize(self):
        if not exists(f'{self.name}_opt.json'):
            print(f'....optimizing {self.name}')
            raise NotImplementedError()

    def get_energy(self):
        print(f'....getting energy of {self.name}')
        raise NotImplementedError()

    def compare_UHF_values(self):
        print(f'....comparing UHF of {self.name}')
        raise NotImplementedError()

    def analyze_cage_geometry(self):
        print(f'....analyzing {self.name}')
        raise NotImplementedError()

    def __str__(self):
        return (
            f'{self.__class__.__name__}'
            f'(name={self.name}, topology={self.topology})'
        )

    def __repr__(self):
        return str(self)


class HetPrism:
    """
    Class that builds and analyses stk.ConstructuedMolecules.

    Represents heteroleptic prism cages and all necessary homoleptic
    cages.

    """
