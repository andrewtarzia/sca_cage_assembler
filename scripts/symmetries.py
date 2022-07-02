#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Classes of symmetries of cages.

Author: Andrew Tarzia

Date Created: 03 Mar 2020

"""


class Symmetry:
    """
    Define a symmetry for a cage topology_graph.

    """

    def __str__(self):
        return (
            f'{self.__class__.__name__}'
            f'{self.definition}'
        )

    def __repr__(self):
        return str(self)


class M8L6_Symmetry(Symmetry):
    """
    Define a symmetry for a M8L6 cage topology_graph.

    """

    def __init__(
        self,
        D_complex,
        L_complex,
        linker,
    ):
        self.D_complex = D_complex
        self.L_complex = L_complex
        self.linker = linker
        self.n_metals = 8
        self.no_vertices = 14

    def d2(self):
        """
        D2 or O symmetry cage.

        Only one distinct symmetry required.

        Delta symmetry complexes at all nodes.

        """

        # All default orientation. All Delta.
        vertex_alignments = {
            i: 0 for i in range(self.no_vertices)
        }
        building_blocks = {
            self.D_complex: range(self.n_metals),
            self.linker: range(self.n_metals, self.no_vertices)
        }

        return {
            'building_blocks': building_blocks,
            'vertex_alignments': vertex_alignments,
            'ratio': (8, 0)
        }

    def th1(self):
        """
        Th symmetry cage 1.

        Two distinct symmetries required.

        """

        # With default rotation pattern.
        # Set metal complex symmetries.
        vertex_alignments = {
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 0, 9: 1, 10: 1, 11: 1, 12: 0, 13: 1
        }
        building_blocks = {
            self.D_complex: (0, 2, 5, 7),
            self.L_complex: (1, 3, 4, 6),
            self.linker: range(self.n_metals, self.no_vertices)
        }

        return {
            'building_blocks': building_blocks,
            'vertex_alignments': vertex_alignments,
            'ratio': (4, 4)
        }

    def th2(self):
        """
        Th symmetry cage 1.

        Two distinct symmetries required.

        """

        # With opposite rotation pattern.
        # Set metal complex symmetries.
        vertex_alignments = {
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 1, 9: 0, 10: 0, 11: 0, 12: 1, 13: 0
        }
        building_blocks = {
            self.D_complex: (0, 2, 5, 7),
            self.L_complex: (1, 3, 4, 6),
            self.linker: range(self.n_metals, self.no_vertices)
        }

        return {
            'building_blocks': building_blocks,
            'vertex_alignments': vertex_alignments,
            'ratio': (4, 4)
        }

    def td(self):
        """
        T symmetry cage.

        One distinct symmetry required.

        All Delta symmetry complexes.

        """

        # With default rotation pattern.
        # All Delta, orientation pattern 1.
        vertex_alignments = {
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 0, 9: 1, 10: 1, 11: 1, 12: 0, 13: 1
        }
        building_blocks = {
            self.D_complex: range(self.n_metals),
            self.linker: range(self.n_metals, self.no_vertices)
        }

        return {
            'building_blocks': building_blocks,
            'vertex_alignments': vertex_alignments,
            'ratio': (8, 0)
        }

    def tl(self):
        """
        T symmetry cage.

        One distinct symmetry required.

        All Lambda symmetry complexes.

        """

        # With default rotation pattern.
        # All Delta, orientation pattern 1.
        vertex_alignments = {
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 0, 9: 1, 10: 1, 11: 1, 12: 0, 13: 1
        }
        building_blocks = {
            self.L_complex: range(self.n_metals),
            self.linker: range(self.n_metals, self.no_vertices)
        }

        return {
            'building_blocks': building_blocks,
            'vertex_alignments': vertex_alignments,
            'ratio': (0, 8)
        }

    def s61(self):
        """
        S_6 symmetry cage 1.

        Two distinct symmetries required.

        Distinct structures have the same complex symmetry patterns
        with different linker rotations.

        """

        # With default rotation pattern.
        # Set metal complex symmetries.
        vertex_alignments = {
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 0, 9: 1, 10: 1, 11: 1, 12: 0, 13: 1
        }
        building_blocks = {
            self.D_complex: (0, 1, 2, 5),
            self.L_complex: (3, 4, 6, 7),
            self.linker: range(self.n_metals, self.no_vertices)
        }

        return {
            'building_blocks': building_blocks,
            'vertex_alignments': vertex_alignments,
            'ratio': (4, 4)
        }

    def s62(self):
        """
        S_6 symmetry cage 2.

        Two distinct symmetries required.

        Distinct structures have the same complex symmetry patterns
        with different linker rotations.

        """

        # With opposite rotation pattern.
        # Set metal complex symmetries.
        vertex_alignments = {
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 1, 9: 0, 10: 0, 11: 0, 12: 1, 13: 0
        }
        building_blocks = {
            self.D_complex: (0, 1, 2, 5),
            self.L_complex: (3, 4, 6, 7),
            self.linker: range(self.n_metals, self.no_vertices)
        }

        return {
            'building_blocks': building_blocks,
            'vertex_alignments': vertex_alignments,
            'ratio': (4, 4)
        }

    def d31(self):
        """
        D_3 symmetry cage 1.

        Two distinct symmetries required.

        Distinct structures have the same complex symmetry patterns
        with different linker rotations.

        Main metal complex symmetry is Delta.

        """

        # With default rotation pattern.
        # Set metal complex symmetries.
        vertex_alignments = {
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 0, 9: 1, 10: 0, 11: 0, 12: 1, 13: 1
        }
        building_blocks = {
            self.D_complex: (0, 2, 3, 4, 5, 6),
            self.L_complex: (1, 7),
            self.linker: range(self.n_metals, self.no_vertices)
        }

        return {
            'building_blocks': building_blocks,
            'vertex_alignments': vertex_alignments,
            'ratio': (6, 2)
        }

    def d32(self):
        """
        D_3 symmetry cage 2.

        Two distinct symmetries required.

        Distinct structures have the same complex symmetry patterns
        with different linker rotations.

        Main metal complex symmetry is Delta.

        """

        # With opposite rotation pattern.
        # Set metal complex symmetries.
        vertex_alignments = {
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 1, 9: 0, 10: 1, 11: 1, 12: 0, 13: 0
        }
        building_blocks = {
            self.D_complex: (0, 2, 3, 4, 5, 6),
            self.L_complex: (1, 7),
            self.linker: range(self.n_metals, self.no_vertices)
        }

        return {
            'building_blocks': building_blocks,
            'vertex_alignments': vertex_alignments,
            'ratio': (6, 2)
        }

    def d31n(self):
        """
        New, homochiral (all delta), D_3 symmetry cage 1.

        Two distinct symmetries required.

        """

        # With default rotation pattern.
        # Set metal complex symmetries.
        vertex_alignments = {
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 0, 9: 1, 10: 0, 11: 0, 12: 1, 13: 1
        }
        building_blocks = {
            self.D_complex: (0, 1, 2, 3, 4, 5, 6, 7),
            self.linker: range(self.n_metals, self.no_vertices)
        }

        return {
            'building_blocks': building_blocks,
            'vertex_alignments': vertex_alignments,
            'ratio': (8, 0)
        }

    def d32n(self):
        """
        New, homochiral (all lambda), D_3 symmetry cage 2.

        Two distinct symmetries required.

        """

        # With opposite rotation pattern.
        # Set metal complex symmetries.
        vertex_alignments = {
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 0, 9: 1, 10: 0, 11: 0, 12: 1, 13: 1
        }
        building_blocks = {
            self.L_complex: (0, 1, 2, 3, 4, 5, 6, 7),
            self.linker: range(self.n_metals, self.no_vertices)
        }

        return {
            'building_blocks': building_blocks,
            'vertex_alignments': vertex_alignments,
            'ratio': (0, 8)
        }

    def s41(self):
        """
        S_4 symmetry cage.

        """

        # With opposite rotation pattern.
        # Set metal complex symmetries.
        vertex_alignments = {
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 1, 9: 1, 10: 0, 11: 1, 12: 0, 13: 0
        }
        building_blocks = {
            self.D_complex: (0, 1, 2, 3),
            self.L_complex: (4, 5, 6, 7),
            self.linker: range(self.n_metals, self.no_vertices)
        }

        return {
            'building_blocks': building_blocks,
            'vertex_alignments': vertex_alignments,
            'ratio': (4, 4)
        }

    def s42(self):
        """
        S_4 symmetry cage.

        """

        # With opposite rotation pattern.
        # Set metal complex symmetries.
        vertex_alignments = {
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 1, 9: 0, 10: 1, 11: 0, 12: 0, 13: 1
        }
        building_blocks = {
            self.D_complex: (0, 1, 2, 3),
            self.L_complex: (4, 5, 6, 7),
            self.linker: range(self.n_metals, self.no_vertices)
        }

        return {
            'building_blocks': building_blocks,
            'vertex_alignments': vertex_alignments,
            'ratio': (4, 4)
        }

    def c2v(self):
        """
        C_2V symmetry cage.

        Only one distinct symmetry required.

        Split metal complex symmetries.

        """

        # With default rotation pattern.
        # Set metal complex symmetries.
        vertex_alignments = {
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 0, 9: 1, 10: 1, 11: 1, 12: 0, 13: 1
        }
        building_blocks = {
            self.D_complex: (1, 2, 4, 7),
            self.L_complex: (0, 3, 5, 6),
            self.linker: range(self.n_metals, self.no_vertices)
        }

        return {
            'building_blocks': building_blocks,
            'vertex_alignments': vertex_alignments,
            'ratio': (4, 4)
        }

    def c2h(self):
        """
        C_2H symmetry cage.

        Only one distinct symmetry required.

        Split metal complex symmetries.

        """

        # With default rotation pattern.
        # Set metal complex symmetries.
        vertex_alignments = {
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 0, 9: 1, 10: 1, 11: 1, 12: 0, 13: 1
        }
        building_blocks = {
            self.D_complex: (0, 1, 2, 3),
            self.L_complex: (4, 5, 6, 7),
            self.linker: range(self.n_metals, self.no_vertices)
        }

        return {
            'building_blocks': building_blocks,
            'vertex_alignments': vertex_alignments,
            'ratio': (4, 4)
        }
