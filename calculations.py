#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for comparing structures.

Author: Andrew Tarzia

Date Created: 15 Mar 2019
"""

import MDAnalysis as mda
from MDAnalysis.analysis.rms import rmsd
import sys


def calculate_RMSD(results_part, structure_name, init_structure_dir):
    '''Calculate RMSD of GFN output structure compared to initial structure.

    Code from James Pegg.

    '''
    # read in initial structure
    initial_structure_file = init_structure_dir+structure_name+'.xyz'
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
