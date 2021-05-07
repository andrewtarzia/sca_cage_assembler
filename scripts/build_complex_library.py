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
import json

import stk
import stko

from molecule_building import (
    custom_fg_factories,
    available_topologies,
    metal_FFs,
    optimize_SCA_complex,
)
from utilities import read_lib


def get_spin_state_energies(complex, name, dict):
    """
    Calculate total electronic energy of each spin state of complex.

    Energy obtained from GFN2-xTB.

    """

    spin_energy_file = f'{name}_spin_energies.json'

    # Iterate over number of unpaired electrons possible.
    spin_energies = {}
    for upe in dict['unpaired_e']:
        xtb_energy = stko.XTBEnergy(
            xtb_path='/home/atarzia/anaconda3/envs/sca_cages/bin/xtb',
            output_dir=f'{name}_{upe}_xtbey',
            charge=dict['total_charge'],
            num_unpaired_electrons=upe,
            unlimited_memory=True,
        )
        energy = xtb_energy.get_energy(complex)
        spin_energies[upe] = energy

    with open(spin_energy_file, 'w') as f:
        json.dump(spin_energies, f)


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
            custom_fg_factories(i)
            for i in ['CNBr_metal', 'CNC_metal']
        ]
        coord_species = stk.BuildingBlock.init_from_file(
            join(ligand_directory, comp['coord_species']),
            functional_groups=ligand_fg_factories
        )

        topology_builder = available_topologies(
            comp['topology']
        )

        complex = topology_builder(metal=metal, ligand=coord_species)
        complex.write(f'{name}.mol')

        # Not interested in unpaired electron checks at this stage.
        # So just select first one.
        range_of_unpaired_e = comp['unpaired_e']
        comp['unpaired_e'] = comp['unpaired_e'][0]

        # Define metal_FFs to use in optimisation.
        custom_metal_FFs = metal_FFs(CN=6)
        complex = optimize_SCA_complex(
            complex=complex,
            name=name,
            dict=comp,
            metal_FFs=custom_metal_FFs
        )

        complex.write(output)

        comp['unpaired_e'] = range_of_unpaired_e
        get_spin_state_energies(
            complex=complex,
            name=name,
            dict=comp,
        )

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
