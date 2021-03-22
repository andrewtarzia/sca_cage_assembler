#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build ligand library.

Author: Andrew Tarzia

Date Created: 13 Jan 2020
"""

import sys
from os.path import exists
import stk
import stko

from molecule_building import (
    metal_FFs,
    custom_fg_factories,
    available_topologies,
)
from utilities import (
    read_lib,
    get_planar_conformer,
)


def build_organics(ligs):

    for name in ligs:
        if ligs[name]['no_metals'] > 0:
            continue
        input_file = f"{ligs[name]['file_loc']}{name}.mol"
        output_file = f'{name}_opt.mol'
        planar_file = f'{name}_planar.mol'
        if exists(output_file):
            continue
        print(f'...building {name}')

        if exists(input_file):
            mol = stk.BuildingBlock.init_from_file(input_file)
        else:
            smi = ligs[name]['smiles']
            mol = stk.BuildingBlock(smiles=smi)
            # Get a planar conformer and optimize crudely.
            if exists(planar_file):
                mol = mol.with_structure_from_file(planar_file)
            else:
                mol = get_planar_conformer(mol)
                mol.write(planar_file)
        mol.write(output_file)


def optimise_metal_centre(name, charge, complex, metal_FF):

    print(f'.......UFF4MOF optimisation of {name}')
    gulp_opt = stko.GulpUFFOptimizer(
        gulp_path='/home/atarzia/software/gulp-5.1/Src/gulp/gulp',
        metal_FF=metal_FF,
        output_dir=f'{name}_uff1'
    )
    gulp_opt.assign_FF(complex)
    complex = gulp_opt.optimize(mol=complex)
    complex.write(f'{name}_uff1.mol')

    print(f'.......XTB optimisation of {name}')
    xtb_opt = stko.XTB(
        xtb_path='/home/atarzia/software/xtb-6.3.1/bin/xtb',
        output_dir=f'{name}_xtb',
        gfn_version=2,
        num_cores=6,
        opt_level='tight',
        charge=charge,
        num_unpaired_electrons=0,
        max_runs=1,
        calculate_hessian=False,
        unlimited_memory=True
    )
    complex = xtb_opt.optimize(mol=complex)
    complex.write(f'{name}_opt.mol')

    return complex


def build_metal_organics(ligs):

    m4_FFs = metal_FFs(CN=4)
    m6_FFs = metal_FFs(CN=6)

    # Iterate over required library.
    for name in ligs:
        if ligs[name]['no_metals'] == 0:
            continue
        comp = ligs[name]
        opt_name = f'{name}_opt.mol'
        if exists(opt_name):
            continue
        print(f'.......building {name}')

        # Build metal atom.
        metal = stk.BuildingBlock(
            smiles=comp['metal_smiles'],
            functional_groups=(
                stk.SingleAtom(stk.Atom(
                    id=0,
                    charge=2,
                    atomic_number=comp['metal_atom_no'],
                ))
                for i in range(comp['no_metal_fgs'])
            ),
            position_matrix=[[0, 0, 0]],
        )

        mc_topo_fn = available_topologies(
            comp['ctopo']
        )

        # Load in organic BB.
        ligand_fg_factories = [
            custom_fg_factories(i)
            for i in comp['organic_FG']
        ]
        organic_BB = stk.BuildingBlock.init_from_file(
            f"{comp['organic_BB']}_opt.mol",
            functional_groups=ligand_fg_factories
        )
        complex = mc_topo_fn(metal, organic_BB)
        complex.write(f'{name}.mol')

        # Optimise metal centre.
        # Assign metal UFF forcefields.
        mffs = m4_FFs if comp['metal_FF'] == 'm4' else m6_FFs
        complex = optimise_metal_centre(
            name=name,
            charge=comp['net_charge'],
            complex=complex,
            metal_FF=mffs
        )
        # Output.
        complex.write(opt_name)


def main():
    if (not len(sys.argv) == 2):
        print("""
Usage: build_ligand_library.py lib_file

    lib_file : (str)
        File containing ligand information (name, smiles, flag)

    """)
        sys.exit()
    else:
        lib_file = sys.argv[1]

    ligs = read_lib(lib_file)

    # Build and optimise all organic molecules in lib.
    build_organics(ligs)

    # Build and optimise all metal containing ligands.
    build_metal_organics(ligs)


if __name__ == "__main__":
    main()
