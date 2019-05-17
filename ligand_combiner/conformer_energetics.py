#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to analyze the energetics and structures of N conformers of all molecules
in an stk.Population()

Author: Andrew Tarzia

Date Created: 03 May 2019

"""

import sys
import logging
import stk
sys.path.insert(0, '/home/atarzia/thesource/')
import analysis


def main():
    if (not len(sys.argv) == 1):
        print("""
    Usage: confomer_energetics.py
        """)
        sys.exit()
    else:
        pass

    logging.info('loading in populations')
    # load in population
    molecule_pop = stk.Population.load('molecules.pop',
                                       member_init=stk.Molecule.from_dict)
    N_mols = len(molecule_pop)
    logging.info('done')
    logging.info('----------------------------------')

    logging.info(f'running energetic analysis on all conformers of {N_mols} molecules:')
    for molecule in molecule_pop:
        N_conf = molecule.mol.GetNumConformers()
        logging.info(f'doing {molecule.name} with {N_conf} conformers')
        analysis.analyze_conformer_NNdist(stk_mol=molecule,
                                          name=molecule.name)
        analysis.analyze_conformer_angles(stk_mol=molecule,
                                          name=molecule.name)
        analysis.analyze_conformer_energies(stk_mol=molecule,
                                            name=molecule.name)

    logging.info('----------------------------------')
    logging.info('done')
    logging.info('----------------------------------')


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='')
    main()
