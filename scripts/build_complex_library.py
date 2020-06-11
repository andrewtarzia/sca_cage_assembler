#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build complex library.

Author: Andrew Tarzia

Date Created: 27 Jan 2020
"""

import sys
from os.path import exists, join

import stk

import molecule_building
from utilities import read_lib


def build_complexes(complexes, ligand_directory):

    for name in complexes:
        comp = complexes[name]
        print(comp)
        output = f'{name}_opt.mol'
        if exists(output):
            continue
        print(f'doing {name}')

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
            molecule_building.custom_fg_factories(i)
            for i in ['CNBr_metal', 'CNC_metal']
        ]
        coord_species = stk.BuildingBlock.init_from_file(
            join(ligand_directory, comp['coord_species']),
            functional_groups=ligand_fg_factories
        )

        topology_builder = molecule_building.available_topologies(
            comp['topology']
        )

        complex = topology_builder(metal=metal, ligand=coord_species)
        complex.write(f'{name}.mol')

        # Not interested in unpaired electron checks at this stage.
        # So just select first one.
        comp['unpaired_e'] = comp['unpaired_e'][0]

        # Define metal_FFs to use in optimisation.
        custom_metal_FFs = molecule_building.metal_FFs(CN=6)
        complex = molecule_building.optimize_SCA_complex(
            complex=complex,
            name=name,
            dict=comp,
            metal_FFs=custom_metal_FFs
        )

        complex.write(output)

    return


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

    compls = read_lib(lib_file)

    # Build and optimise all organic molecules in lib.
    build_complexes(compls, ligand_directory)


if __name__ == "__main__":
    main()
