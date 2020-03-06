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
from os.path import exists, join
import glob
from rdkit.Chem import AllChem as rdkit
import stk

import cage_building


def read_lib(lib_file):
    """
    Read lib file.

    Returns dictionary.

    """

    with open(lib_file, 'r') as f:
        lib = json.load(f)

    return lib


def build_cages(
    complexes,
    prisms,
    ligand_directory,
    complex_directory
):

    for name in prisms:
        pris = prisms[name]
        compl_names = pris['corners'].strip(')(').split(', ')
        print(compl_names)
        comps = {i: complexes[i] for i in compl_names}
        print(pris)
        print(comps)
        het_cage = cage_building.HetPrism(
            name=name,
            prism_dict=pris,
            complex_dicts=comps,
            ligand_dir=ligand_directory,
            complex_dir=complex_directory
        )
        # print(het_cage)
        print('tob', het_cage.cages_to_build)
        for C in het_cage.cages_to_build:
            print(C)
            C.build()
            default_free_e = C.free_electron_options[0]
            print(C.free_electron_options, default_free_e)
            input()
            C.optimize(free_e=default_free_e)
            sys.exit()
        sys.exit()


def main():
    first_line = (
        'Usage: build_cage_library.py prism_lib_file compl_lib_file'
        'lig_directory compl_directory'
    )
    if (not len(sys.argv) == 5):
        print(f"""
{first_line}

    compl_lib_file : (str)
        File containing complex information (XXXXX).

    prism_lib_file : (str)
        File containing prism information (XXXXX).

    lig_directory : (str)
        Directory with required ligand structures.

    compl_directory : (str)
        Directory with required complex structures.

    """)
        sys.exit()
    else:
        compl_lib_file = sys.argv[1]
        prism_lib_file = sys.argv[2]
        ligand_directory = sys.argv[3]
        compl_directory = sys.argv[4]

    print(f'reading {prism_lib_file}')
    prisms = read_lib(prism_lib_file)
    print(f'reading {compl_lib_file}')
    compls = read_lib(compl_lib_file)

    # Build and optimise all organic molecules in lib.
    build_cages(compls, prisms, ligand_directory, compl_directory)

    sys.exit()


if __name__ == "__main__":
    main()
