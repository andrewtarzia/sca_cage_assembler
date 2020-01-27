#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build complex library.

Author: Andrew Tarzia

Date Created: 27 Jan 2020
"""
import sys
import json
def read_complex_lib(lib_file):
    """
    Read complex lib file.

    Returns dictionary of format:

    ligs[name] = (smiles, flag)

    """
    with open(lib_file, 'r') as f:
        compls = json.load(f)

    print(compls)
    return compls

def main():
    if (not len(sys.argv) == 3):
        print("""
Usage: build_complex_library.py lib_file ligand_directory

    lib_file : (str)
        File containing complex information (XXXXX).

    ligand_directory : (str)
        Directory with required ligand structures.

    """)
        sys.exit()
    else:
        lib_file = sys.argv[1]
        ligand_directory = sys.argv[2]


if __name__ == "__main__":
    main()
