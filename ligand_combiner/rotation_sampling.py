#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build N conformers of stk.Polymers() and determine GFN energy as fn
of a rotation coordination.

Author: Andrew Tarzia

Date Created: 03 May 2019

"""

import sys
from numpy.linalg import norm
from numpy import asarray, cos, radians
from itertools import product
from pandas import read_csv
from rdkit.Chem import AllChem as Chem
from os.path import join, isfile
from os import remove
from glob import glob
from stk import Population, Molecule
from Combiner import get_molecule, get_geometrical_properties, Combination
sys.path.insert(0, '/home/atarzia/thesource/')
from stk_functions import build_population
from rdkit_functions import mol_list2grid
from analysis import plot_all_pair_info, output_analysis


def main():
    if (not len(sys.argv) == 8):
        print("""
    Usage: molecule_builing.py subset rebuild N bond_mean bond_std energy_tol angle_tol
        subset (str) - 'all' to build all possible molecules,
            'clever' to build molecules from bloch2017
        rebuild (str) - 't' if you want to rebuild all molecules, 'f' to load populations
        N (int) - number of conformers to use
        bond_mean (float) - mean value of bond distance to use in candidate selection
        bond_std (float) - std deviation value of bond distance to use in candidate selection
        energy_tol (float) - max kJ/mol over min energy conformer to allow
        angle_tol (float) - tolerance to use on angle matching
        """)
        sys.exit()
    else:
        subset = sys.argv[1]
        rebuild = sys.argv[2]
        N = int(sys.argv[3])
        bond_mean = float(sys.argv[4])
        bond_std = float(sys.argv[5])
        energy_tol = float(sys.argv[6])
        angle_tol = float(sys.argv[7])

    # set molecules to load

    # iterate over molecule

    # load in stk polymer

    # build N conformers

    # determine geometrical coordinate

    # run GFN SPE -- get total energy

    # plot energy vs coordinate


if __name__ == "__main__":
    main()
