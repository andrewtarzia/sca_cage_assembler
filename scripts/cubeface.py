#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module of classes for face construction.

Author: Andrew Tarzia

Date Created: 25 Jan 2021

"""

from itertools import combinations
from scipy.spatial.distance import euclidean
import numpy as np

import stk
from stk.utilities import (
    get_acute_vector,
    get_plane_normal,
)


from atools import get_atom_distance


class _FaceVertex(stk.Vertex):
    """
    Represents a vertex of a :class:`.CubeFace`.

    """

    def __init__(
        self,
        id,
        position,
        use_neighbor_placement=True,
        aligner_edge=0,
    ):
        """
        Initialize a :class:`._CageVertex`.

        Parameters
        ----------
        id : :class:`int`
            The id of the vertex.

        position : :class:`tuple` of :class:`float`
            The position of the vertex.

        use_neighbor_placement : :class:`bool`, optional
            If ``True``, the position of the vertex will be updated
            based on the neighboring functional groups.

        aligner_edge : :class:`int`, optional
            The edge which is used to align the :class:`.BuildingBlock`
            placed on the vertex. The first :class:`.FunctionalGroup`
            is rotated such that it lies exactly on this
            :class:`.Edge`. Must be between ``0`` and the number of
            edges the vertex is connected to.

        """

        self._use_neighbor_placement = use_neighbor_placement
        self._aligner_edge = aligner_edge
        super().__init__(id, position)

    def clone(self):
        clone = super().clone()
        clone._aligner_edge = self._aligner_edge
        clone._use_neighbor_placement = self._use_neighbor_placement
        return clone

    def _with_aligner_edge(self, aligner_edge):
        """
        Modify the instance.

        """

        self._aligner_edge = aligner_edge
        return self

    def with_aligner_edge(self, aligner_edge):
        """
        Return a clone with a different `aligner_edge`.

        Parameters
        ----------
        aligner_edge : :class:`int`
            The aligner edge of the clone.

        Returns
        -------
        :class:`._CageVertex`
            The clone. Has the same type as the original instance.

        """

        return self.clone()._with_aligner_edge(aligner_edge)

    def use_neighbor_placement(self):
        """
        ``True`` if the position should be updated based on neighbors.

        Returns
        -------
        :class:`bool`
            ``True`` if the position of the vertex should be updated
            based on the positions of functional groups on neighboring
            vertices.

        """

        return self._use_neighbor_placement

    @classmethod
    def init_at_center(cls, id, vertices):
        """
        Initialize a :class:`._CageVertex` in the middle of `vertices`.

        Parameters
        ----------
        id : :class:`int`
            The id of the initialized vertex.

        vertices : :class:`tuple` of :class:`.Vertex`
            The vertices at whose center this one needs to be.

        Returns
        -------
        :class:`._CageVertex`
            The new vertex.

        """

        return cls(
            id=id,
            position=(
                sum(vertex.get_position() for vertex in vertices)
                / len(vertices)
            ),
        )

    def get_aligner_edge(self):
        """
        Return the aligner edge of the vertex.

        Returns
        -------
        :class:`int`
            The aligner edge.

        """

        return self._aligner_edge

    def __str__(self):
        return (
            f'Vertex(id={self._id}, '
            f'position={self._position.tolist()}, '
            f'aligner_edge={self._aligner_edge})'
        )


class _MetalVertex(_FaceVertex):
    def place_building_block(self, building_block, edges):
        assert (
            building_block.get_num_functional_groups() == 1
        ), (
            f'{building_block} needs to have exactly 1 functional '
            'groups but has '
            f'{building_block.get_num_functional_groups()}.'
        )
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )

        # Get normal to placer plane (NCCN).
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        normal = building_block.get_plane_normal(
            atom_ids=building_block.get_placer_ids(),
        )
        normal = get_acute_vector(
            reference=core_centroid - self._position,
            vector=normal,
        )

        # Align it along Z direction. [CHECK THIS DIRECTION]
        building_block = building_block.with_rotation_between_vectors(
            start=normal,
            target=[0, 0, 1],
            origin=self._position,
        )

        # Align fg vector with edge.
        fg, = building_block.get_functional_groups()
        fg_start_centroid = building_block.get_centroid(
            atom_ids=[i for i in fg.get_placer_ids()][:2],
        )
        fg_end_centroid = building_block.get_centroid(
            atom_ids=[i for i in fg.get_placer_ids()][2:],
        )
        edge_centroid = (
            sum(edge.get_position() for edge in edges) / len(edges)
        )
        return building_block.with_rotation_between_vectors(
            start=fg_end_centroid - fg_start_centroid,
            target=edge_centroid - self._position,
            origin=self._position,
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):

        return {
            fg_id: edge.get_id() for fg_id, edge in enumerate(edges)
        }


class _LinkerVertex(_FaceVertex):

    def _get_building_block_long_axis(self, building_block):

        for fg in building_block.get_functional_groups():
            if len(list(fg.get_bonder_ids())) > 1:
                raise ValueError(
                    f'{building_block} has functional groups with more'
                    ' than 1 binder.'
                )

        binder_atom_ids = [
            list(fg.get_bonder_ids())
            for fg in building_block.get_functional_groups()
        ]
        binder_atom_dists = sorted(
            [
                (idx1, idx2, get_atom_distance(
                    building_block,
                    idx1,
                    idx2
                ))
                for idx1, idx2 in combinations(binder_atom_ids, r=2)
            ],
            key=lambda a: a[2]
        )

        # Can assume the ordering of the binder atom distances:
        # 0, 1: short vectors
        # 2, 3: long vectors
        # 4, 5: diagonal vectors
        # This fails when the molecule is not sufficiently anisotropic,
        # at which point it will not matter.
        short_vector_fg_1 = (
            binder_atom_dists[0][0], binder_atom_dists[0][1]
        )
        short_vector_fg_2 = (
            (binder_atom_dists[1][0], binder_atom_dists[1][1])
            if (
                binder_atom_dists[1][0] not in short_vector_fg_1 and
                binder_atom_dists[1][1] not in short_vector_fg_1
            ) else
            (binder_atom_dists[2][0], binder_atom_dists[2][1])
        )

        short_pair_1_centroid = building_block.get_centroid(
            atom_ids=(i[0] for i in short_vector_fg_1),
        )
        short_pair_2_centroid = building_block.get_centroid(
            atom_ids=(i[0] for i in short_vector_fg_2),
        )

        long_axis = short_pair_2_centroid - short_pair_1_centroid
        return long_axis

    def place_building_block(self, building_block, edges):
        assert (
            building_block.get_num_functional_groups() == 4
        ), (
            f'{building_block} needs to have 4 functional '
            'groups but has '
            f'{building_block.get_num_functional_groups()}.'
        )
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )

        edge_centroid = (
            sum(edge.get_position() for edge in edges) / len(edges)
        )
        edge_normal = get_acute_vector(
            reference=edge_centroid,
            vector=get_plane_normal(
                points=np.array([
                    edge.get_position() for edge in edges
                ]),
            ),
        )

        fg_bonder_centroid = building_block.get_centroid(
            atom_ids=next(
                building_block.get_functional_groups()
            ).get_placer_ids(),
        )
        edge_position = edges[self._aligner_edge].get_position()
        building_block = (
            building_block.with_rotation_to_minimize_angle(
                start=fg_bonder_centroid - self._position,
                target=edge_position - edge_centroid,
                axis=edge_normal,
                origin=self._position,
            )
        )

        # Flatten wrt to xy plane.
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        normal = building_block.get_plane_normal(
            atom_ids=building_block.get_placer_ids(),
        )
        normal = get_acute_vector(
            reference=core_centroid - self._position,
            vector=normal,
        )
        building_block = building_block.with_rotation_between_vectors(
            start=normal,
            target=[0, 0, 1],
            origin=self._position,
        )

        # Align long axis of molecule (defined by FG centroid) with
        # X axis.
        long_axis_vector = self._get_building_block_long_axis(
            building_block
        )
        building_block = (
            building_block.with_rotation_to_minimize_angle(
                start=long_axis_vector,
                target=[1, 0, 0],
                axis=edge_normal,
                origin=self._position,
            )
        )
        return building_block.get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):

        def fg_distance(edge):
            return euclidean(edge.get_position(), fg_position)

        # For each FG, get the closest edge.
        mapping = {}
        for fg_id, fg in enumerate(
            building_block.get_functional_groups()
        ):
            fg_position = building_block.get_centroid(
                fg.get_placer_ids()
            )
            edges = sorted(edges, key=fg_distance)
            mapping[fg_id] = edges[0].get_id()

        return mapping


class CubeFace(stk.cage.Cage):
    """
    Represents a cage topology graph.

    Metal building blocks with three functional groups are
    required for this topology.

    Ligand building blocks with four functional groups are required for
    this topology.

    When using a :class:`dict` for initialization, a
    :class:`.BuildingBlock` needs to be assigned to each of the
    following numbers:

        | metals (3 functional groups): 0 to 3
        | ligands (4 functional groups): 4

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        _MetalVertex(
            id=0,
            position=[1, 1, 0],
            use_neighbor_placement=False,
        ),
        _MetalVertex(
            id=1,
            position=[1, -1, 0],
            use_neighbor_placement=False,
        ),
        _MetalVertex(
            id=2,
            position=[-1, -1, 0],
            use_neighbor_placement=False,
        ),
        _MetalVertex(
            id=3,
            position=[-1, 1, 0],
            use_neighbor_placement=False,
        ),
    )

    _vertex_prototypes = (
        *_vertex_prototypes,

        _LinkerVertex(
            id=4,
            position=[0, 0, 0],
            use_neighbor_placement=False,
        ),
    )

    _edge_prototypes = (
        stk.Edge(
            id=0,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[4],
            position=[0.5, 0.5, 0],
        ),
        stk.Edge(
            id=1,
            vertex1=_vertex_prototypes[1],
            vertex2=_vertex_prototypes[4],
            position=[0.5, -0.5, 0],
        ),
        stk.Edge(
            id=2,
            vertex1=_vertex_prototypes[2],
            vertex2=_vertex_prototypes[4],
            position=[-0.5, -0.5, 0],
        ),
        stk.Edge(
            id=3,
            vertex1=_vertex_prototypes[3],
            vertex2=_vertex_prototypes[4],
            position=[-0.5, 0.5, 0],
        ),
    )
