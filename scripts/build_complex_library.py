#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build complex library.

Author: Andrew Tarzia

"""

import sys
import os
import stk
import logging

from molecule_building import (
    custom_fg_factories,
    available_topologies,
    metal_FFs,
    optimize_SCA_complex,
)
from env_set import read_envset_json, read_json

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)


def build_complexes(complexes, environment):

    for name in complexes:
        comp = complexes[name]
        output = environment["complex_dir"] / f'{name}_opt.mol'
        if os.path.exists(output):
            continue
        logging.info(f'doing {name}')

        # Build metal atom.
        metal = stk.BuildingBlock(
            smiles=comp['metal_smiles'],
            functional_groups=(
                stk.SingleAtom(stk.Atom(
                    id=0,
                    charge=2,
                    atomic_number=comp['metal_atom_no'],
                ))
                for i in range(6)
            ),
            position_matrix=[[0, 0, 0]],
        )

        ligand_fg_factories = [
            custom_fg_factories(i)
            for i in ['CNBr_metal', 'CNC_metal']
        ]
        coord_species = stk.BuildingBlock.init_from_file(
            path=str(environment["ligand_dir"] / comp['coord_species']),
            functional_groups=ligand_fg_factories,
        )

        topology_builder = available_topologies(comp['topology'])

        compl_bb = topology_builder(metal=metal, ligand=coord_species)
        compl_bb.write(environment["complex_dir"] / f'{name}.mol')

        # Not interested in unpaired electron checks at this stage.
        # So just select first one.
        range_of_unpaired_e = comp['unpaired_e']
        comp['unpaired_e'] = comp['unpaired_e'][0]

        # Define metal_FFs to use in optimisation.
        custom_metal_FFs = metal_FFs(CN=6)
        compl_bb = optimize_SCA_complex(
            compl=compl_bb,
            name=name,
            dict=comp,
            metal_FFs=custom_metal_FFs,
            environment=environment,
        )

        compl_bb.write(output)

        comp['unpaired_e'] = range_of_unpaired_e


def main():
    if (not len(sys.argv) == 2):
        print("""
Usage: build_complex_library.py env_set_file

    env_set_file : (str)
        File containing environment set file.


    """)
        sys.exit()
    else:
        env_set_file = sys.argv[1]

    environment = read_envset_json(env_set_file)
    compls = read_json(environment["complex_file"])

    # Build and optimise all organic molecules in lib.
    build_complexes(compls, environment)


if __name__ == "__main__":
    main()
