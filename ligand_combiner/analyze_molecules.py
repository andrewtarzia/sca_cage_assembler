#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to analyze all built molecules from DB of core, ligands and linkers.

Author: Andrew Tarzia

Date Created: 17 Apr 2019

"""

import sys
import pickle
import numpy as np
import pandas as pd
from rdkit.Chem import AllChem as Chem
from os.path import join
import matplotlib.pyplot as plt
import stk
from stk import Population
from Combiner import get_molecule, get_geometrical_properties, Combination
sys.path.insert(0, '/home/atarzia/thesource/')
from stk_functions import build_population
from rdkit_functions import mol_list2grid


def main():
    if (not len(sys.argv) == 3):
        print("""
    Usage: analyze_molecules.py pair_data angle_tol
        pair_data (str) - pickle file with pair data
        angle_tol (float) - tolerance to use on angle matching
        """)
        sys.exit()
    else:
        pair_data = sys.argv[1]
        angle_tol = sys.argv[2]

    proj_dir = '/home/atarzia/projects/ligand_combiner/'
    core_dir = proj_dir + 'cores/'
    liga_dir = proj_dir + 'ligands/'
    link_dir = proj_dir + 'linkers/'
    mole_dir = proj_dir + 'molecules/'

    # load in pair data
    with open(join(mole_dir, pair_data), 'rb') as f:
        all_pairs = pickle.load(f)
    print(len(all_pairs))
