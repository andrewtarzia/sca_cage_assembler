#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build ligand library.

Author: Andrew Tarzia

"""

import sys
from os.path import exists
import stk
import bbprep
import logging

from env_set import read_envset_json, read_json

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)


def build_organics(ligs, environment):

    for name in ligs:
        if ligs[name]["no_metals"] > 0:
            continue

        planar_file = environment["ligand_dir"] / f"{name}_planar.mol"
        if exists(planar_file):
            continue
        logging.info(f"...building {name}")

        smi = ligs[name]["smiles"]
        mol = stk.BuildingBlock(smiles=smi)
        # Get a planar conformer.
        if not exists(planar_file):
            ensemble = bbprep.generators.ETKDG(50).generate_conformers(mol)
            process = bbprep.Planarfy(
                ensemble=ensemble,
                selector=bbprep.selectors.AllNonHSelector(),
            )
            min_molecule = process.get_minimum()
            min_molecule.molecule.write(planar_file)


def main():
    if not len(sys.argv) == 2:
        print(
            """
Usage: build_ligand_library.py env_set_file

    env_set_file : (str)
        File containing environment set file.

    """
        )
        sys.exit()
    else:
        env_set_file = sys.argv[1]

    environment = read_envset_json(env_set_file)

    ligs = read_json(environment["lib_file"])

    # Build and optimise all organic molecules in lib.
    build_organics(ligs, environment)


if __name__ == "__main__":
    main()
