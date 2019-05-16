#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to iterate through all mol files in DIR and get UFF energy.

Author: Andrew Tarzia

Date Created: 16 May 2019

"""

import sys
import logging
import stk
import glob
from rdkit.Chem import AllChem as Chem
sys.path.insert(0, '/home/atarzia/thesource/')
import Combiner
import stk_f


def main():
    if (not len(sys.argv) == 1):
        print("""
    Usage: get_struct_energy.py
        """)
        sys.exit()
    else:
        pass

    molecules = sorted(glob.glob('*.mol'))
    for mol in molecules:
        stk_mol = stk_f.load_StructUnitX(file=mol, X=0)
        energy = Combiner.get_energy(stk_mol=stk_mol, conformer=-1, FF='UFF')
        logging.info(f'molecule {mol} has UFF energy = {energy} kJ/mol')


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='')
    main()
