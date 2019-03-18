#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to get RMSD between two structures.

Author: Andrew Tarzia

Date Created: 08 Mar 2019

"""

import MDAnalysis as mda
from MDAnalysis.analysis.rms import rmsd
import sys
sys.path.insert(0, '/home/atarzia/thesource/')


if __name__ == "__main__":
    if (not len(sys.argv) == 3):
        print("""
Usage: get_RMSD.py struct1 struct2
    struct1: template structure
    struct2: structure to compare
""")
        sys.exit()
    else:
        xyz1 = sys.argv[1]
        xyz2 = sys.argv[2]

    # read in template structure
    initial_structure_file = xyz1
    ref = mda.Universe(initial_structure_file)
    # read in new structure
    new_structure_file = xyz2
    mobile = mda.Universe(new_structure_file)
    RMSD = rmsd(mobile.atoms.positions, ref.atoms.positions, center=True,
                superposition=True)
    print(xyz1, xyz2, '--->', RMSD)
