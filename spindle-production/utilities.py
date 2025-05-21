"""Utilities module."""

import json

import numpy as np
from scipy.spatial.distance import euclidean
from stk import get_acute_vector, get_plane_normal


def read_lib(lib_file):
    """Read lib file."""
    with open(lib_file, "rb") as f:
        return json.load(f)


def get_atom_distance(molecule, atom1_id, atom2_id):
    """Return the distance between atom1 and atom2.

    Parameters
    ----------
    molecule : :class:`stk.Molecule`

    atom1_id : :class:`int`
        The id of atom1.

    atom2_id : :class:`int`
        The id of atom2.

    Returns:
    -------
    :class:`float`
        The euclidean distance between two atoms.

    """
    position_matrix = molecule.get_position_matrix()

    distance = euclidean(
        u=position_matrix[atom1_id], v=position_matrix[atom2_id]
    )

    return float(distance)


def reorient_linker(molecule):
    target_coords = (
        np.array([1, 1, 0]),
        np.array([1, -1, 0]),
        np.array([-1, -1, 0]),
        np.array([-1, 1, 0]),
    )
    centroid_pos = np.array([0, 0, 0])
    molecule = molecule.with_centroid(
        position=centroid_pos,
        atom_ids=molecule.get_placer_ids(),
    )

    edge_centroid = sum(target_coords) / len(target_coords)
    edge_normal = get_acute_vector(
        reference=edge_centroid,
        vector=get_plane_normal(
            points=np.array(target_coords),
        ),
    )

    fg_bonder_centroid = molecule.get_centroid(
        atom_ids=next(molecule.get_functional_groups()).get_placer_ids(),
    )
    edge_position = target_coords[0]
    molecule = molecule.with_rotation_to_minimize_angle(
        start=fg_bonder_centroid - centroid_pos,
        target=edge_position - edge_centroid,
        axis=edge_normal,
        origin=centroid_pos,
    )

    # Flatten wrt to xy plane.
    core_centroid = molecule.get_centroid(
        atom_ids=molecule.get_core_atom_ids(),
    )
    normal = molecule.get_plane_normal(
        atom_ids=molecule.get_placer_ids(),
    )
    normal = get_acute_vector(
        reference=core_centroid - centroid_pos,
        vector=normal,
    )
    molecule = molecule.with_rotation_between_vectors(
        start=normal,
        target=[0, 0, 1],
        origin=centroid_pos,
    )

    # Align long axis of molecule (defined by deleter atoms) with
    # y axis.
    long_axis_vector = molecule.get_long_axis()
    edge_centroid = sum(target_coords) / len(target_coords)
    edge_normal = get_acute_vector(
        reference=edge_centroid,
        vector=get_plane_normal(
            points=np.array(target_coords),
        ),
    )
    molecule = molecule.with_rotation_to_minimize_angle(
        start=long_axis_vector,
        target=[1, 0, 0],
        axis=edge_normal,
        origin=centroid_pos,
    )
    return molecule
