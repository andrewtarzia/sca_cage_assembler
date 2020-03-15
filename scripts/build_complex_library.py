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

import stk

import molecule_building


def read_complex_lib(lib_file):
    """
    Read complex lib file.

    Returns dictionary.

    """

    with open(lib_file, 'r') as f:
        compls = json.load(f)

    return compls


def build_complexes(complexes, ligand_directory):

    for name in complexes:
        comp = complexes[name]
        print(comp)
        output = f'{name}_opt.mol'
        jsonoutput = f'{name}_opt.json'
        if exists(jsonoutput):
            continue
        print(f'doing {name}')

        # Build metal atom.
        metal = molecule_building.build_metal(
            metal_smiles=comp['metal_smiles'],
            no_fgs=6
        )
        # Define binding atom and binding FG.
        binding_atom = molecule_building.build_atom(
            'N',
            FG='metal_bound_N'
        )
        binding_fgs = ['metal_bound_N']
        # Build initial metal centre for all complexes.
        # Always a six coordinate, octahedral complex with N atoms.
        metal_centre = molecule_building.build_metal_centre(
            metal=metal,
            topology=stk.metal_centre.Octahedral(),
            binding_atom=binding_atom,
            return_FG=binding_fgs
        )
        metal_centre.write(f'{name}_metal_centre.mol')

        coord_species = stk.BuildingBlock.init_from_file(
            join(ligand_directory, comp['coord_species']),
            ['CNC_metal', 'CNBr_metal']
        )
        coord_species = molecule_building.order_FGs(
            mol=coord_species,
            order=['CNBr_metal', 'CNC_metal']
        )

        topology = molecule_building.available_topologies(
            comp['topology']
        )

        complex = molecule_building.build_SCA_complex(
            metal_centre=metal_centre,
            complex_top=topology,
            bidentate_ligand=coord_species
        )
        complex.write(f'{name}.mol')
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
        complex.dump(jsonoutput)

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

    print(f'reading {lib_file}')
    compls = read_complex_lib(lib_file)

    # Build and optimise all organic molecules in lib.
    build_complexes(compls, ligand_directory)


if __name__ == "__main__":
    main()
