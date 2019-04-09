#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for comparing structures.

Author: Andrew Tarzia

Date Created: 15 Mar 2019
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.rms import rmsd
import sys


def calculate_RMSD(results_part, structure_name, init_structure_dir):
    '''Calculate RMSD of GFN output structure compared to initial structure.

    Code from James Pegg.

    '''
    # read in initial structure
    initial_structure_file = init_structure_dir + structure_name + '.xyz'
    ref = mda.Universe(initial_structure_file)
    # read in new structure
    new_structure_file = 'xtbopt.xyz'
    mobile = mda.Universe(new_structure_file)
    # RMSD of all atoms
    RMSD = rmsd(mobile.atoms.positions,
                ref.atoms.positions,
                center=True, superposition=True)
    results_part['RMSD'] = RMSD
    return results_part


def calc_formation_energy(prod, react):
    '''Calculate formation energy of 'A' in a.u. from 2 lists of energies.

    Reaction formation energy == sum(product energy) - sum(reactant energy)

    Keyword arguments:
        prod (list) - list of product energies
        react (list) - list of reactant energies

    Returns:
        RFE (float) - Reaction formation energy in a.u.

    '''
    RFE = sum(prod) - sum(react)
    return RFE


def get_formation_energies(data, ff='OPLS'):
    '''Calculate formation energies based on prod - react energies.

    '''
    # by definition
    if ff == 'OPLS':
        H2O_energy = 0.0
    else:
        print('need water energy - need to implement check for this in data')
        sys.exit('exitting')
    form_energies = []
    stoich = {'1p1': {'bb1': 1, 'bb2': 1, 'h2o': 3},
              '4p4': {'bb1': 4, 'bb2': 4, 'h2o': 12},
              '2p3': {'bb1': 2, 'bb2': 3, 'h2o': 6},
              '4p6': {'bb1': 4, 'bb2': 6, 'h2o': 12},
              '4p62': {'bb1': 4, 'bb2': 6, 'h2o': 12},
              '6p9': {'bb1': 6, 'bb2': 9, 'h2o': 18}}
    for i, row in data.iterrows():
        FE = (row.cage_ey + stoich[row.topo]['h2o'] * H2O_energy) - (row.bb1_ey * stoich[row.topo]['bb1'] + row.bb2_ey * stoich[row.topo]['bb2'])
        form_energies.append(FE)
    return form_energies


def unit_vector(vector):
    """ Returns the unit vector of the vector.

    https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python/13849249#13849249
    """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793

    https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python/13849249#13849249
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def get_dihedral(pt1, pt2, pt3, pt4):
    '''Calculate the dihedral (-pi to pi) between four points using Praxeolitic formula
    1 sqrt, 1 cross product

    From: https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    (new_dihedral(p))
    '''
    p0 = pt1
    p1 = pt2
    p2 = pt3
    p3 = pt4

    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))
