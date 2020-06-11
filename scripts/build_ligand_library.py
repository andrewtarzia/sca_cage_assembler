#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build ligand library.

Author: Andrew Tarzia

Date Created: 13 Jan 2020
"""

import sys
from os import mkdir
from os.path import exists
from rdkit.Chem import AllChem as rdkit
import stk
import stko

import atools
import molecule_building
from utilities import read_lib, calculate_binding_AR


def build_organics(ligs):

    optimizer = stko.XTB(
        xtb_path='/home/atarzia/software/xtb-6.3.0/bin/xtb',
        output_dir='lig_opt_',
        gfn_version=2,
        num_cores=6,
        opt_level='extreme',
        max_runs=1,
        calculate_hessian=False,
        unlimited_memory=True
    )

    for name in ligs:
        if ligs[name]['no_metals'] > 0:
            continue
        input_file = f"{ligs[name]['file_loc']}{name}.mol"
        output_file = f'{name}_opt.mol'
        if exists(output_file):
            continue
        print(f'...building {name}')

        if exists(input_file):
            mol = stk.BuildingBlock.init_from_file(input_file)
        else:
            smi = ligs[name]['smiles']
            mol = stk.BuildingBlock(smiles=smi)
        optimizer.optimize(mol)
        mol.write(output_file)


def optimise_metal_centre(name, charge, complex, metal_FF):

    print(f'.......UFF4MOF optimisation of {name}')
    gulp_opt = stko.GulpMetalOptimizer(
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

    m4_FFs = molecule_building.metal_FFs(CN=4)
    m6_FFs = molecule_building.metal_FFs(CN=6)

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

        mc_topo_fn = molecule_building.available_topologies(
            comp['ctopo']
        )

        # Load in organic BB.
        ligand_fg_factories = [
            molecule_building.custom_fg_factories(i)
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


def output_2d_images(ligs):

    if not exists('built_ligands/'):
        mkdir('built_ligands')

    # Draw 2D representation of all built molecules.
    mols = []
    names = []
    for name in ligs:
        comp = ligs[name]
        opt_name = f'{name}_opt.mol'
        if ligs[name]['no_metals'] == 0:
            MOL = rdkit.MolFromMolFile(opt_name)
            MOL.RemoveAllConformers()
            mols.append(MOL)
            AR = calculate_binding_AR(
                stk.BuildingBlock.init_from_file(
                    opt_name,
                    functional_groups=[stk.BromoFactory()]
                )
            )
            label = f'{name}'
            if AR is not None:
                label = f'{label}: {round(AR,2)}'
            names.append(label)
        else:
            MOL = rdkit.MolFromMolFile(opt_name)
            MOL.RemoveAllConformers()
            mols.append(MOL)
            top_str = comp['ctopo']
            label = f'{name}\n{top_str}'
            AR = calculate_binding_AR(
                stk.BuildingBlock.init_from_file(
                    opt_name,
                    functional_groups=[stk.BromoFactory()]
                )
            )
            if AR is not None:
                label = f'{label}: {round(AR,2)}'
            names.append(label)

    # Draw 2D representation of all built molecules.
    atools.mol_list2grid(
        molecules=mols,
        names=names,
        filename='built_ligands/built_ligands',
        mol_per_row=3,
        maxrows=3,
        subImgSize=(250, 200)
    )

    return


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

    # Produce image of all built molecules.
    output_2d_images(ligs)


if __name__ == "__main__":
    main()
