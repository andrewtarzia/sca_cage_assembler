#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Classes of face sets for defined topologies and symmetries.

Author: Andrew Tarzia

Date Created: 05 Nov 2020

"""


class FaceSets:
    """
    Define face sets for a cage topology_graph.

    """

    def __str__(self):
        return (f'{self.__class__.__name__}')

    def __repr__(self):
        return str(self)


class M8L6_FaceSets(FaceSets):
    """
    Define a symmetry for a M8L6 cage topology_graph.

    """

    def __init__(self, symmetry_string):
        self.symmetry_string = symmetry_string
        self.sets = ['001', '100', '010', '00-1', '-100', '0-10']
        self.vertices = self._define_vertices(symmetry_string)
        self.connected = self._define_connected(symmetry_string)
        self.opposite = self._define_opposite()

    def _define_opposite(self):
        return {
            '001': '00-1',
            '100': '-100',
            '010': '0-10',
            '00-1': '001',
            '-100': '100',
            '0-10': '010',
        }

    def _define_vertices(self, symmetry_string):
        if symmetry_string in ['o1', 't1']:
            return {
                '001': (0, 1, 2, 3),
                '100': (0, 1, 4, 5),
                '010': (0, 3, 4, 7),
                '00-1': (4, 5, 6, 7),
                '-100': (2, 3, 6, 7),
                '0-10': (1, 2, 5, 6),
            }
        elif symmetry_string in ['th1', 'th2']:
            return {
                '001': (0, 1, 4, 5),
                '100': (0, 2, 4, 6),
                '010': (0, 3, 5, 6),
                '00-1': (2, 3, 6, 7),
                '-100': (1, 3, 5, 7),
                '0-10': (1, 2, 4, 7),
            }
        elif symmetry_string in ['s61', 's62']:
            return {
                '001': (0, 1, 2, 4),
                '100': (0, 1, 3, 5),
                '010': (0, 4, 5, 7),
                '00-1': (3, 5, 6, 7),
                '-100': (2, 4, 6, 7),
                '0-10': (1, 2, 3, 6),
            }
        elif symmetry_string in ['d31', 'd32']:
            return {
                '001': (0, 1, 2, 6),
                '100': (0, 3, 4, 6),
                '010': (0, 2, 3, 7),
                '00-1': (3, 4, 5, 7),
                '-100': (2, 1, 5, 7),
                '0-10': (1, 4, 5, 6),
            }
        elif symmetry_string in ['c2v']:
            return {
                '001': (0, 1, 4, 5),
                '100': (0, 2, 4, 6),
                '010': (2, 3, 4, 5),
                '00-1': (2, 3, 6, 7),
                '-100': (1, 3, 5, 7),
                '0-10': (0, 1, 6, 7),
            }
        elif symmetry_string in ['c2h']:
            return {
                '001': (0, 1, 2, 3),
                '100': (0, 1, 4, 5),
                '010': (0, 3, 4, 7),
                '00-1': (4, 5, 6, 7),
                '-100': (2, 3, 6, 7),
                '0-10': (1, 2, 5, 6),
            }

    def _define_connected(self, symmetry_string):
        if symmetry_string in ['o1', 't1']:
            return {
                '001': ((0, 1), (1, 2), (2, 3), (3, 0)),
                '100': ((0, 1), (1, 5), (5, 4), (4, 0)),
                '010': ((0, 3), (3, 7), (7, 4), (4, 0)),
                '00-1': ((4, 5), (5, 6), (6, 7), (7, 4)),
                '-100': ((2, 3), (3, 7), (7, 6), (6, 2)),
                '0-10': ((1, 2), (2, 6), (6, 5), (5, 1)),
            }
        elif symmetry_string in ['th1', 'th2']:
            return {
                '001': ((0, 4), (4, 1), (1, 5), (5, 0)),
                '100': ((0, 4), (4, 2), (2, 6), (6, 0)),
                '010': ((0, 6), (6, 3), (3, 5), (5, 0)),
                '00-1': ((3, 6), (6, 2), (2, 7), (7, 3)),
                '-100': ((1, 5), (5, 3), (3, 7), (7, 1)),
                '0-10': ((1, 4), (4, 2), (2, 7), (7, 1)),
            }
        elif symmetry_string in ['s61', 's62']:
            return {
                '001': ((0, 1), (1, 2), (2, 4), (4, 0)),
                '100': ((0, 5), (5, 3), (3, 1), (1, 0)),
                '010': ((0, 4), (4, 7), (7, 5), (5, 0)),
                '00-1': ((3, 6), (6, 7), (7, 5), (5, 3)),
                '-100': ((4, 2), (2, 6), (6, 7), (7, 4)),
                '0-10': ((1, 2), (2, 6), (6, 3), (3, 1)),
            }
        elif symmetry_string in ['d31', 'd32']:
            return {
                '001': ((0, 6), (6, 1), (1, 2), (2, 0)),
                '100': ((0, 3), (3, 4), (4, 6), (6, 0)),
                '010': ((0, 2), (2, 7), (7, 3), (3, 0)),
                '00-1': ((7, 3), (3, 4), (4, 5), (5, 7)),
                '-100': ((2, 1), (1, 5), (5, 7), (7, 2)),
                '0-10': ((6, 1), (1, 5), (5, 4), (4, 6)),
            }
        elif symmetry_string in ['c2v']:
            return {
                '001': ((0, 1), (1, 5), (5, 4), (4, 0)),
                '100': ((0, 6), (6, 2), (2, 4), (4, 0)),
                '010': ((2, 3), (3, 5), (5, 4), (4, 2)),
                '00-1': ((2, 3), (3, 7), (7, 6), (6, 2)),
                '-100': ((1, 7), (7, 3), (3, 5), (5, 1)),
                '0-10': ((0, 1), (1, 7), (7, 6), (6, 0)),
            }
        elif symmetry_string in ['c2h']:
            return {
                '001': ((0, 1), (1, 2), (2, 3), (3, 0)),
                '100': ((0, 1), (1, 5), (5, 4), (4, 0)),
                '010': ((0, 3), (3, 7), (7, 4), (4, 0)),
                '00-1': ((4, 5), (5, 6), (6, 7), (7, 4)),
                '-100': ((2, 3), (3, 7), (7, 6), (6, 2)),
                '0-10': ((1, 2), (2, 6), (6, 5), (5, 1)),
            }