#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build and analyse faces from ligand library.

Author: Andrew Tarzia

Date Created: 19 Oct 2020

"""

import sys
from os.path import exists, join
from glob import glob
import numpy as np

import stk
from scipy.spatial.distance import euclidean

from stk.molecular.topology_graphs.utilities import (
    _FunctionalGroupSorter,
    _EdgeSorter,
)
from stk.utilities import (
    get_acute_vector,
    get_plane_normal,
    normalize_vector,
)
import stko

from molecule_building import metal_FFs


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
    def place_building_block(self, building_block, edges):
        assert (
            building_block.get_num_functional_groups() > 2
        ), (
            f'{building_block} needs to have more than 2 functional '
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
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        placer_centroid = building_block.get_centroid(
            atom_ids=building_block.get_placer_ids(),
        )
        building_block = building_block.with_rotation_between_vectors(
            start=get_acute_vector(
                reference=core_centroid - placer_centroid,
                vector=building_block.get_plane_normal(
                    atom_ids=building_block.get_placer_ids(),
                ),
            ),
            target=edge_normal,
            origin=self._position,
        )
        fg_bonder_centroid = building_block.get_centroid(
            atom_ids=next(
                building_block.get_functional_groups()
            ).get_placer_ids(),
        )
        edge_position = edges[self._aligner_edge].get_position()
        return building_block.with_rotation_to_minimize_angle(
            start=fg_bonder_centroid - self._position,
            target=edge_position - edge_centroid,
            axis=edge_normal,
            origin=self._position,
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        # The idea is to order the functional groups in building_block
        # by their angle with the vector running from the placer
        # centroid to fg0, going in the clockwise direction.
        # The edges are also ordered by their angle with the vector
        # running from the edge centroid to the aligner_edge,
        # going in the clockwise direction.
        #
        # Once the fgs and edges are ordered, zip and assign them.

        fg_sorter = _FunctionalGroupSorter(building_block)
        edge_sorter = _EdgeSorter(
            edges=edges,
            aligner_edge=edges[self._aligner_edge],
            axis=fg_sorter.get_axis(),
        )
        return {
            fg_id: edge.get_id()
            for fg_id, edge in zip(
                fg_sorter.get_items(),
                edge_sorter.get_items(),
            )
        }


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


def load_complex(filename):

    fgfactory = stk.SmartsFunctionalGroupFactory(
        smarts='[#7X3]~[#6]~[#6]~[#7X3]~[#35]',
        bonders=(3, ),
        deleters=(4, ),
        placers=(0, 1, 2, 3),
    )

    name = filename.replace('.mol', '')
    # Need to define more than one placer id for complexes to ensure
    # alignment -- use the NCCN plane (attached to the Br) to define
    # the orientation of the complex.
    complex = stk.BuildingBlock.init_from_file(
        filename,
        functional_groups=[fgfactory]
    )

    return name, complex


def optimize_complex(complex, name):

    opt_name = f'{name}_opt.mol'
    if exists(opt_name):
        return complex.with_structure_from_file(opt_name)
    else:
        print(f'doing UFF4MOF optimisation for {name}')
        gulp_opt = stko.GulpUFFOptimizer(
            gulp_path='/home/atarzia/software/gulp-5.1/Src/gulp/gulp',
            metal_FF=metal_FFs(CN=6),
            output_dir=f'{name}_uff1'
        )
        gulp_opt.assign_FF(complex)
        complex = gulp_opt.optimize(mol=complex)
        complex.write(f'{name}_uff1.mol')

        print(f'doing xTB optimisation for {name}')
        xtb_opt = stko.XTB(
            xtb_path='/home/atarzia/software/xtb-6.3.1/bin/xtb',
            output_dir=f'{name}_xtb',
            gfn_version=2,
            num_cores=6,
            opt_level='tight',
            charge=2,
            num_unpaired_electrons=0,
            max_runs=1,
            calculate_hessian=False,
            unlimited_memory=True
        )
        complex = xtb_opt.optimize(mol=complex)
        complex.write(f'{name}_opt.mol')


def load_ligands(directory):

    ligands = {}
    for lig in glob(join(directory, f'quad2*_opt.mol')):
        l_name = lig.replace(directory, '').replace('_opt.mol', '')
        ligands[l_name] = stk.BuildingBlock.init_from_file(
            lig,
            functional_groups=[stk.BromoFactory()],
        )

    return ligands


def build_face(
    face_name,
    lig_structure,
    del_complex,
    lam_complex,
    face_topo
):

    face_file = f'{face_name}.mol'

    face = stk.ConstructedMolecule(
        topology_graph=CubeFace(
            building_blocks={
                del_complex: face_topo['d_pos'],
                lam_complex: face_topo['l_pos'],
                lig_structure: (4, ),
            },
            vertex_alignments=face_topo['va'],
        )
    )

    face.write(face_file)

def get_planar_conformer(molecule):
    cids, confs = build_conformers(
        mol=molecule,
        N=100,
        ETKDG_version='v3'
    )
    print(f'getting optimal conformer...')
    min_plane_dev = 100000000
    min_cid = -10

    new_molecule = molecule.clone()

    for cid in cids:

        # Update stk_mol to conformer geometry.
        new_molecule = update_from_rdkit_conf(
            stk_mol=new_molecule,
            rdk_mol=confs,
            conf_id=cid
        )

        plane_dev = calculate_molecule_planarity(new_molecule)
        if plane_dev < min_plane_dev:
            min_cid = cid
            min_plane_dev = plane_dev
            molecule = update_from_rdkit_conf(
                stk_mol=molecule,
                rdk_mol=confs,
                conf_id=min_cid
            )

    return molecule


def planarfy(ligands):
    """
    Get the most planar conformer of each ligand.

    This is done by determining the ETKDG conformer with the smallest
    plane deviation from its plane of best fit.

    """

    new_ligands = {}

    for ligand in ligands:
        planar_file = f'{ligand}_planar.mol'
        if exists(planar_file):
            opt_lig = ligands[ligand].with_structure_from_file(
                planar_file
            )
        else:
            print(f'doing {ligand}...')
            opt_lig = get_planar_conformer(ligands[ligand])
            opt_lig.write(planar_file)
        new_ligands[ligand] = opt_lig

    return new_ligands


def check_ligand_orientation(mol):
    """
    Check the ligand orientation assuming four functional groups.

    The orientation here assumes the ligand makes a pseudo rectangle
    with its FGs, and that we want the ``short axis'' to be defined by
    FG(0)-FG(1) vector.

    This requires each FG to only have one bonder atom.

    """

    if mol.get_num_functional_groups() != 4:
        raise ValueError(
            f'{mol} has {mol.get_num_functional_groups()}, '
            'but needs 4'
        )

    print(mol)
    mol.write('temp.xyz')
    needs_rotation = False

    for fg in mol.get_functional_groups():
        if len(list(fg.get_bonder_ids())) > 1:
            raise ValueError(
                f'{mol} has functional groups with more than 1 binder.'
            )

    binder_atom_ids = [
        list(fg.get_bonder_ids())
        for fg in mol.get_functional_groups()
    ]
    binder_atom_dists = sorted(
        [
            (idx1, idx2, get_atom_distance(mol, idx1, idx2))
            for idx1, idx2 in combinations(binder_atom_ids, r=2)
        ],
        key=lambda a: a[2]
    )
    print(binder_atom_ids)
    print(binder_atom_dists)

    # Can assume the ordering of the binder atom distances:
    # 0, 1: short vectors
    # 2, 3: long vectors
    # 4, 5: diagonal vectors
    # This fails when the molecule is not sufficiently anisotropic, at
    # which point it will not matter.
    short_vector_fgs = [
        (binder_atom_dists[0][0], binder_atom_dists[0][1]),
        (binder_atom_dists[1][0], binder_atom_dists[1][1]),
    ]
    long_vector_fgs = [
        (binder_atom_dists[2][0], binder_atom_dists[2][1]),
        (binder_atom_dists[3][0], binder_atom_dists[3][1]),
    ]
    print(short_vector_fgs, long_vector_fgs)

    chk1 = binder_atom_ids[0] in short_vector_fgs[0]
    chk2 = binder_atom_ids[1] in short_vector_fgs[0]
    if (chk1 and chk2):
        print('short is right')
        needs_rotation = False
    else:
        input('wrong?')
        needs_rotation = True

    return needs_rotation


def main():
    first_line = (
        'Usage: face_analysis.py '
        'lig_directory'
    )
    if (not len(sys.argv) == 2):
        print(f"""
{first_line}

    ligand_directory : (str)
        Directory with required ligand structures.

    """)
        sys.exit()
    else:
        ligand_directory = sys.argv[1]

    # Define two metal building blocks (lambda, delta).
    del_name, del_complex = load_complex('cl1_zn_oct_del_face.mol')
    lam_name, lam_complex = load_complex('cl1_zn_oct_lam_face.mol')

    # Optimise both complexes.
    del_complex = optimize_complex(del_complex, del_name)
    lam_complex = optimize_complex(lam_complex, lam_name)

    # Load in each ligand structure.
    ligands = load_ligands(ligand_directory)

    # Get most planar conformer of ligand.
    ligands = planarfy(ligands)

    face_topologies = {
        '1': {
            'va': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
            'ratio': (2, 2),
            'd_pos': (1, 3),
            'l_pos': (0, 2),
        },
        '2': {
            'va': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
            'ratio': (2, 2),
            'd_pos': (0, 2),
            'l_pos': (1, 3),
        },
        '3': {
            'va': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
            'ratio': (4, 0),
            'd_pos': (),
            'l_pos': (0, 1, 2, 3),
        },
        '4': {
            'va': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
            'ratio': (3, 1),
            'd_pos': (3, ),
            'l_pos': (0, 1, 2),
        },
        '5': {
            'va': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
            'ratio': (3, 1),
            'd_pos': (2, ),
            'l_pos': (0, 1, 3, ),
        },
    }

    # Build and optimise five face options per ligand.
    face_properties = {}
    for lig in sorted(ligands):
        print(f'doing {lig}...')
        lig_structure = ligands[lig]
        # Check if ligand needs rotation.
        # needs_rotation = check_ligand_orientation(lig_structure)
        # print(needs_rotation)
        # continue
        # needs_rotation = False
        face_properties[lig] = {}
        # Build each face topology.
        for face_t in face_topologies:
            final_topology_dict = face_topologies[face_t]
            # print(face_topologies[face_t], needs_rotation)
            # if needs_rotation:
            #     final_topology_dict = face_topologies[face_t]
            #     final_topology_dict['va'][4: 1]
            # else:
            #     final_topology_dict = face_topologies[face_t]
            # print(final_topology_dict)
            # sys.exit()
            face_name = f'{lig}_{face_t}'
            build_face(
                face_name=face_name,
                lig_structure=lig_structure,
                del_complex=del_complex,
                lam_complex=lam_complex,
                face_topo=face_topologies[face_t],
            )
            continue
            optimize_face(face_name)

            # Measure properties.
            get_face_properties(face_name)
            face_properties[lig][face_name] = []

    sys.exit()


if __name__ == '__main__':
    main()
