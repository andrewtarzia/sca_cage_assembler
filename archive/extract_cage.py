#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to extract cage molecule from CIF.

Author: Andrew Tarzia

Date Created: 15 Mar 2019

"""

import sys


def main():
    if (not len(sys.argv) == 2):
        print("""
Usage: extract_cage.py CIF
    CIF: file (.cif) to analyze and add pseudo atoms to.
    """)
        sys.exit()
    else:
        file = sys.argv[1]

    file_prefix = file.replace('.cif', '')
    pdb_file, ASE_structure = convert_CIF_2_PDB(file)
    if pdb_file is None and ASE_structure is None:
        sys.exit()

    # Rebuild system.
    rebuilt_structure = modularize(file=pdb_file)
    if rebuilt_structure is None:
        # Handle pyWindow failure.
        sys.exit(f'pyWindow failure on {pdb_file}')

    for molecule in rebuilt_structure.molecules:
        print('Analysing molecule {0} out of {1}'.format(
            molecule + 1, len(rebuilt_structure.molecules))
        )
        mol = rebuilt_structure.molecules[molecule]
        # skip very small structures -- likely ions
        if mol.no_of_atoms < 20:
            continue

        print(mol.full_analysis(), '\n')
        # Each molecule can be saved separately
        pdb_output = file_prefix + "_{0}.pdb".format(molecule)
        mol.dump_molecule(
            pdb_output,
            include_coms=False,
            override=True)
        _, _ = convert_PDB_2_XYZ(pdb_output)


if __name__ == "__main__":
    main()
