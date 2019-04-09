#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to draw the precusors files in a directory using RDKit.

Author: Andrew Tarzia

Date Created: 09 Apr 2019
"""

from glob import glob
import sys
from rdkit.Chem import AllChem as Chem
sys.path.insert(0, '/home/atarzia/thesource/')
from rdkit_functions import mol_list2grid


def main():
    """Run script.

    """
    if (not len(sys.argv) == 3):
        print("""
Usage: draw_precursors.py suffix output
    suffix (str) - file type to read structures from
    output (str) - name of image file
    """)
        sys.exit()
    else:
        suffix = sys.argv[1]
        output = sys.argv[2]

    molecules = []
    for file in glob('*' + suffix):
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromMolFile(file)))
        molecules.append(mol)

    mol_list2grid(molecules=molecules, names=None, filename=output,
                  mol_per_row=3, subImgSize=(200, 200))


if __name__ == "__main__":
    main()
