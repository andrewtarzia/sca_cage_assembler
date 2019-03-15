#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for comparing structures.

Author: Andrew Tarzia

Date Created: 15 Mar 2019
"""

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd


def calculate_RMSD(results_part, structure_name, init_structure_dir):
    '''Calculate RMSD of GFN output structure compared to initial structure.

    Code from James Pegg.

    '''
    # read in initial structure
    initial_structure_file = init_structure_dir+file+'.xyz'
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
