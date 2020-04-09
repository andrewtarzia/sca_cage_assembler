#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Modules defining and building the Cage class and subclasses.

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
        n_metals,
        topo,
        no_vertices
    ):
        self.D_complex = D_complex
        self.L_complex = L_complex
        self.linker = linker
        self.n_metals = n_metals
        self.topo = topo
        self.no_vertices = no_vertices

    def o1(self):
        """
        O symmetry cage.

        Only one distinct symmetry required.

        Delta symmetry complexes at all nodes.

        """

        # All default orientation.
        orient_1 = self.topo(vertex_alignments={
            i: 0 for i in range(self.no_vertices)
        })
        linker_verts_1 = orient_1.vertices[self.n_metals:]
        # All default orientation. All Delta.
        return (
            orient_1,
            {
                self.D_complex: orient_1.vertices[:self.n_metals],
                self.linker: linker_verts_1
            },
            (8, 0)
        )

    def th1(self):
        """
        Th symmetry cage 1.

        Two distinct symmetries required.

        """

        # With default rotation pattern.
        th_orient_1 = self.topo(vertex_alignments={
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 0, 9: 1, 10: 1, 11: 1, 12: 0, 13: 1
        })

        # Set metal complex symmetries.
        iter = [0, 1, 0, 1, 1, 0, 1, 0]
        th_D_verts_1 = [
            v for i, v in enumerate(
                th_orient_1.vertices[:self.n_metals]
            )
            if iter[i] == 0
        ]
        th_L_verts_1 = [
            v for i, v in enumerate(
                th_orient_1.vertices[:self.n_metals]
            )
            if iter[i] == 1
        ]

        return (
            th_orient_1,
            {
                self.D_complex: th_D_verts_1,
                self.L_complex: th_L_verts_1,
                self.linker: th_orient_1.vertices[self.n_metals:]
            },
            (4, 4)
        )

    def th2(self):
        """
        Th symmetry cage 1.

        Two distinct symmetries required.

        """

        # With opposite rotation pattern.
        th_orient_2 = self.topo(vertex_alignments={
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 1, 9: 0, 10: 0, 11: 0, 12: 1, 13: 0
        })

        iter = [0, 1, 0, 1, 1, 0, 1, 0]
        th_D_verts_2 = [
            v for i, v in enumerate(
                th_orient_2.vertices[:self.n_metals]
            )
            if iter[i] == 0
        ]
        th_L_verts_2 = [
            v for i, v in enumerate(
                th_orient_2.vertices[:self.n_metals]
            )
            if iter[i] == 1
        ]

        # With rotated pattern.
        return (
            th_orient_2,
            {
                self.D_complex: th_D_verts_2,
                self.L_complex: th_L_verts_2,
                self.linker: th_orient_2.vertices[self.n_metals:]
            },
            (4, 4)
        )

    def t1(self):
        """
        T symmetry cage 1.

        One distinct symmetry required.

        All Delta symmetry complexes.

        """

        # With default rotation pattern.
        th_orient_1 = self.topo(vertex_alignments={
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 0, 9: 1, 10: 1, 11: 1, 12: 0, 13: 1
        })

        # All Delta, orientation pattern 1.
        return (
            th_orient_1,
            {
                self.D_complex: th_orient_1.vertices[:self.n_metals],
                self.linker: th_orient_1.vertices[self.n_metals:]
            },
            (8, 0)
        )

    def s61(self):
        """
        S_6 symmetry cage 1.

        Two distinct symmetries required.

        Distinct structures have the same complex symmetry patterns
        with different linker rotations.

        """

        # With default rotation pattern.
        th_orient_1 = self.topo(vertex_alignments={
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 0, 9: 1, 10: 1, 11: 1, 12: 0, 13: 1
        })

        iter = [0, 0, 0, 1, 1, 0, 1, 1]
        s6_D_verts = [
            v for i, v in enumerate(
                th_orient_1.vertices[:self.n_metals]
            )
            if iter[i] == 0
        ]
        s6_L_verts = [
            v for i, v in enumerate(
                th_orient_1.vertices[:self.n_metals]
            )
            if iter[i] == 1
        ]

        # Complex orientation 1, orientation pattern 1.
        return (
            th_orient_1,
            {
                self.D_complex: s6_D_verts,
                self.L_complex: s6_L_verts,
                self.linker: th_orient_1.vertices[self.n_metals:]
            },
            (4, 4)
        )

    def s62(self):
        """
        S_6 symmetry cage 2.

        Two distinct symmetries required.

        Distinct structures have the same complex symmetry patterns
        with different linker rotations.

        """

        # With opposite rotation pattern.
        th_orient_2 = self.topo(vertex_alignments={
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 1, 9: 0, 10: 0, 11: 0, 12: 1, 13: 0
        })

        iter = [0, 0, 0, 1, 1, 0, 1, 1]
        s6_D_verts = [
            v for i, v in enumerate(
                th_orient_2.vertices[:self.n_metals]
            )
            if iter[i] == 0
        ]
        s6_L_verts = [
            v for i, v in enumerate(
                th_orient_2.vertices[:self.n_metals]
            )
            if iter[i] == 1
        ]

        # Complex orientation 1, orientation pattern 2.
        return (
            th_orient_2,
            {
                self.D_complex: s6_D_verts,
                self.L_complex: s6_L_verts,
                self.linker: th_orient_2.vertices[self.n_metals:]
            },
            (4, 4)
        )

    def d31(self):
        """
        D_3 symmetry cage 1.

        Two distinct symmetries required.

        Distinct structures have the same complex symmetry patterns
        with different linker rotations.

        Main metal complex symmetry is Delta.

        """

        # With default rotation pattern.
        d3_orient_1 = self.topo(vertex_alignments={
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 0, 9: 1, 10: 0, 11: 0, 12: 1, 13: 1
        })

        d3_iter_1 = [0, 1, 0, 0, 0, 0, 0, 1]

        d3_D_verts_1 = [
            v for i, v in enumerate(
                d3_orient_1.vertices[:self.n_metals]
            )
            if d3_iter_1[i] == 0
        ]
        d3_L_verts_1 = [
            v for i, v in enumerate(
                d3_orient_1.vertices[:self.n_metals]
            )
            if d3_iter_1[i] == 1
        ]
        return (
            d3_orient_1,
            {
                self.D_complex: d3_D_verts_1,
                self.L_complex: d3_L_verts_1,
                self.linker: d3_orient_1.vertices[self.n_metals:]
            },
            (6, 2)
        )

    def d32(self):
        """
        D_3 symmetry cage 2.

        Two distinct symmetries required.

        Distinct structures have the same complex symmetry patterns
        with different linker rotations.

        Main metal complex symmetry is Delta.

        """

        # With opposite rotation pattern.
        d3_orient_2 = self.topo(vertex_alignments={
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 1, 9: 0, 10: 1, 11: 1, 12: 0, 13: 0
        })

        d3_iter_1 = [0, 1, 0, 0, 0, 0, 0, 1]

        d3_D_verts_1 = [
            v for i, v in enumerate(
                d3_orient_2.vertices[:self.n_metals]
            )
            if d3_iter_1[i] == 0
        ]
        d3_L_verts_1 = [
            v for i, v in enumerate(
                d3_orient_2.vertices[:self.n_metals]
            )
            if d3_iter_1[i] == 1
        ]
        return (
            d3_orient_2,
            {
                self.D_complex: d3_D_verts_1,
                self.L_complex: d3_L_verts_1,
                self.linker: d3_orient_2.vertices[self.n_metals:]
            },
            (6, 2)
        )

    def c2v(self):
        """
        C_2V symmetry cage.

        Only one distinct symmetry required.

        Split metal complex symmetries.

        """

        # With default rotation pattern.
        orient = self.topo(vertex_alignments={
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 0, 9: 1, 10: 1, 11: 1, 12: 0, 13: 1
        })

        # Set metal complex symmetries.
        iter = [1, 0, 0, 1, 0, 1, 1, 0]
        D_verts = [
            v for i, v in enumerate(
                orient.vertices[:self.n_metals]
            )
            if iter[i] == 0
        ]
        L_verts = [
            v for i, v in enumerate(
                orient.vertices[:self.n_metals]
            )
            if iter[i] == 1
        ]

        return (
            orient,
            {
                self.D_complex: D_verts,
                self.L_complex: L_verts,
                self.linker: orient.vertices[self.n_metals:]
            },
            (4, 4)
        )

    def c2h(self):
        """
        C_2H symmetry cage.

        Only one distinct symmetry required.

        Split metal complex symmetries.

        """

        # With default rotation pattern.
        orient = self.topo(vertex_alignments={
            # Metals.
            0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0,
            # Linkers.
            8: 0, 9: 1, 10: 1, 11: 1, 12: 0, 13: 1
        })

        # Set metal complex symmetries.
        iter = [1, 1, 1, 1, 0, 0, 0, 0]
        D_verts = [
            v for i, v in enumerate(
                orient.vertices[:self.n_metals]
            )
            if iter[i] == 0
        ]
        L_verts = [
            v for i, v in enumerate(
                orient.vertices[:self.n_metals]
            )
            if iter[i] == 1
        ]

        return (
            orient,
            {
                self.D_complex: D_verts,
                self.L_complex: L_verts,
                self.linker: orient.vertices[self.n_metals:]
            },
            (4, 4)
        )
