#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build all molecules from DB of core, ligands and linkers.

Author: Andrew Tarzia

Date Created: 05 May 2019

"""

import sys
from stk import Population, Molecule
sys.path.insert(0, '/home/atarzia/thesource/')
from analysis import output_analysis, get_all_pairs


def main():
    if (not len(sys.argv) == 7):
        print("""
    Usage: analyze_molecules_pair.py bond_mean bond_std energy_tol angle_tol mol1 mol2
        bond_mean (float) - mean value of bond distance to use in candidate selection
        bond_std (float) - std deviation value of bond distance to use in candidate selection
        energy_tol (float) - max kJ/mol over min energy conformer to allow
        angle_tol (float) - tolerance to use on angle matching
        mol1 (str) - prefix of molecule to analyze (ABCBA_XX OR ABA_XX)
        mol2 (str) - prefix of molecule to analyze (ABCBA_XX OR ABA_XX)
        """)
        sys.exit()
    else:
        bond_mean = float(sys.argv[1])
        bond_std = float(sys.argv[2])
        energy_tol = float(sys.argv[3])
        angle_tol = float(sys.argv[4])
        mol1 = sys.argv[5]
        mol2 = sys.argv[6]

    print('loading in populations')
    # load in population
    molecule_pop = Population.load('molecules.pop',
                                   member_init=Molecule.from_dict)
    print('done')
    print('----------------------------------')

    print('obtaining properties for all pairs')
    all_pairs = get_all_pairs(molecule_pop=molecule_pop,
                              settings={'bond_mean': bond_mean,
                                        'bond_std': bond_std,
                                        'energy_tol': energy_tol},
                              mol_pair=(mol1, mol2))

    print('done')
    print(len(all_pairs))
    print('----------------------------------')
    # do analysis
    print('doing all analysis')
    output_analysis(molecule_pop=molecule_pop, pair_data=all_pairs,
                    angle_tol=angle_tol, energy_tol=energy_tol,
                    mol_pair=(mol1, mol2))
    print('done')
    print('----------------------------------')


if __name__ == "__main__":
    main()
