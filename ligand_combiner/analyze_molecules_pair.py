#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build all molecules from DB of core, ligands and linkers.

Author: Andrew Tarzia

Date Created: 05 May 2019

"""

import sys
import logging
import stk
sys.path.insert(0, '/home/atarzia/thesource/')
import analysis


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

    logging.info('loading in populations')
    # load in population
    molecule_pop = stk.Population.load('molecules.pop',
                                       member_init=stk.Molecule.from_dict)
    logging.info('done')
    logging.info('----------------------------------')

    logging.info('obtaining properties for all pairs')
    all_pairs = analysis.get_all_pairs(molecule_pop=molecule_pop,
                                       settings={'bond_mean': bond_mean,
                                                 'bond_std': bond_std,
                                                 'energy_tol': energy_tol},
                                       mol_pair=(mol1, mol2))

    logging.info('----------------------------------')
    logging.info('done')
    logging.info(f'number of pairs: {len(all_pairs)}')
    logging.info('----------------------------------')
    # do analysis
    logging.info('doing all analysis')
    analysis.output_analysis(molecule_pop=molecule_pop, pair_data=all_pairs,
                             angle_tol=angle_tol, energy_tol=energy_tol,
                             mol_pair=(mol1, mol2))
    logging.info('done')
    logging.info('----------------------------------')


if __name__ == "__main__":
    # logging.basicConfig(level=logging.DEBUG, format='%(levelname)s-%(message)s')
    # logging.debug(f'Debug mode!')
    logging.basicConfig(level=logging.INFO, format='')
    main()
