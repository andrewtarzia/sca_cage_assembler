#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Classes of symmetries of cages.

Author: Andrew Tarzia

Date Created: 03 Mar 2020

"""

from itertools import product


def all_m8l6face_symmetries(
    topo,
    D_complex,
    L_complex,
    linker,
    check_orientation,
    no_vertices,
    rotatable_vertices,
    complex_vertices
):
    """
    Define a list of all possible M8L6 symmetries.

    """
    symm_list = {}
    topologies_to_build = {}
    if check_orientation:
        # Need to define a list of orientations to use by
        # iterating over the vertex alignments of each
        # vertex.Here we have assumed that each ligand can
        # take twoorientations (original, and 90deg
        # rotation about face plane).
        iteration = product(
            [0, 1],
            repeat=len(rotatable_vertices)
        )
        for i in iteration:
            v_align = {i: 0 for i in range(no_vertices)}
            for j, rv in enumerate(rotatable_vertices):
                v_align[rv] = i[j]
            v_align_string = ''.join([
                str(i) for i in list(v_align.values())
            ])[len(rotatable_vertices):]
            topologies_to_build[v_align_string] = topo(
                vertex_alignments=v_align
            )
    else:
        # Use only a single orientation topology.
        v_align = {i: 0 for i in range(no_vertices)}
        v_align_string = ''.join(list(v_align.values()))
        topologies_to_build[v_align_string] = topo(
            vertex_alignments=v_align
        )
    print(f'{len(topologies_to_build)} topologies')

    n_metals = 8
    # Iterate over all face orientations and complex
    # symmetries.
    for topo in topologies_to_build:
        topo_f = topologies_to_build[topo]

        # For each ratio, must define all the possible
        # places for each complex symmetry.
        # Need to define a list of bb vertex inputs of the
        # two symmetry complexes that produces all
        # possible placementvariation on the cage metal
        # vertices.
        iteration = product(
            [0, 1], repeat=len(complex_vertices)
        )
        for iter in iteration:
            D_verts = [
                v for i, v in enumerate(
                    topo_f.vertices[:n_metals]
                )
                if iter[i] == 0
            ]
            L_verts = [
                v for i, v in enumerate(
                    topo_f.vertices[:n_metals]
                )
                if iter[i] == 1
            ]
            ratio = (len(D_verts), len(L_verts))
            linker_verts = topo_f.vertices[n_metals:]
            bb_vert = {
                D_complex: D_verts,
                L_complex: L_verts,
                linker: linker_verts
            }
            bb_vert_string = ''.join([
                str(i) for i in iter
            ])
            name_string = bb_vert_string+topo
            symm_list[name_string] = (
                topo_f,
                bb_vert,
                ratio
            )
    print(f'{len(symm_list)} symmetries')
    return symm_list


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

    def t(self):
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


class M4L4_Symmetry(Symmetry):
    """
    Define a symmetry for a M4L4 cage topology_graph.

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
        self.n_metals = 4
        self.no_vertices = 8

    def base(self):
        """
        Cage with default orientations and delta symmetry.

        Only one distinct symmetry required.

        Delta symmetry complexes at all nodes.

        """

        # All default orientation. All Delta.
        vertex_alignments = {
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0,
            # Linkers.
            4: 0, 5: 0, 6: 0, 7: 0,
        }
        building_blocks = {
            self.D_complex: range(self.n_metals),
            self.linker: range(self.n_metals, self.no_vertices)
        }

        return {
            'building_blocks': building_blocks,
            'vertex_alignments': vertex_alignments,
            'ratio': (4, 0)
        }


class M6L2L3_Symmetry(Symmetry):
    """
    Define a symmetry for a M6L2L3Prism cage topology_graph.

    """

    def __init__(
        self,
        D_complex,
        L_complex,
        linker3,
        linker4,
    ):
        self.D_complex = D_complex
        self.L_complex = L_complex
        self.linker3 = linker3
        self.linker4 = linker4
        self.n_metals = 6
        self._tritopic_vertices = 2
        self._tetratopic_vertices = 3
        self.no_vertices = 11

    def base(self):
        """
        Cage with default orientations and delta symmetry.

        Only one distinct symmetry required.

        Delta symmetry complexes at all nodes.

        """

        # All default orientation. All Delta.
        vertex_alignments = {
            i: 0 for i in range(self.no_vertices)
        }
        building_blocks = {
            self.D_complex: range(self.n_metals),
            self.linker3: range(
                self.n_metals,
                self.n_metals+self._tritopic_vertices,
            ),
            self.linker4: range(
                self.n_metals+self._tritopic_vertices,
                self.no_vertices,
            ),
        }

        return {
            'building_blocks': building_blocks,
            'vertex_alignments': vertex_alignments,
            'ratio': (6, 0)
        }
