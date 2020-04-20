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

import atools
import molecule_building
from utilities import read_lib, calculate_binding_AR


def build_organics(ligs):

    optimizer = stk.XTB(
        xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
        output_dir='lig_opt_',
        gfn_version=2,
        num_cores=6,
        opt_level='extreme',
        max_runs=1,
        calculate_hessian=False,
        unlimited_memory=True
    )

    for name in ligs:
        if ligs[name]['is_organic'] is False:
            continue
        input_file = f"{ligs[name]['file_loc']}{name}.mol"
        output_file = f'{name}_opt.mol'
        jsonoutput = f'{name}_opt.json'
        if exists(jsonoutput):
            continue
        print(f'...building {name}')

        if exists(input_file):
            mol = stk.BuildingBlock.init_from_file(input_file)
        else:
            smi = ligs[name]['smiles']
            mol = stk.BuildingBlock(smiles=smi)
        optimizer.optimize(mol)
        mol.write(output_file)
        mol.dump(jsonoutput)


def binding_atoms(string):
    """
    Defines available binding atoms for metal-containing ligands.

    Returns stk.BuildingBlock.

    """

    dictionary = {
        'metal_bound_N': molecule_building.build_atom(
            'N',
            FG='metal_bound_N'
        ),
        'metal_bound_O': molecule_building.build_atom(
            'O',
            FG='metal_bound_O'
        )
    }

    if string not in dictionary.keys():
        raise KeyError(
            f'{string} not in dictionary.\n'
            f'available binding atoms: {dictionary.keys()}'
        )

    return dictionary[string]


def optimise_metal_centre(name, charge, complex, metal_FF):

    print(f'.......UFF4MOF optimisation of {name}')
    gulp_opt = stk.GulpMetalOptimizer(
        gulp_path='/home/atarzia/software/gulp-5.1/Src/gulp/gulp',
        metal_FF=metal_FF,
        output_dir=f'{name}_uff1'
    )
    gulp_opt.assign_FF(complex)
    gulp_opt.optimize(mol=complex)
    complex.write(f'{name}_uff1.mol')
    complex.dump(f'{name}_uff1.json')

    print(f'.......XTB optimisation of {name}')
    xtb_opt = stk.XTB(
        xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
        output_dir=f'{name}',
        gfn_version=2,
        num_cores=6,
        opt_level='tight',
        charge=charge,
        num_unpaired_electrons=0,
        max_runs=1,
        calculate_hessian=False,
        unlimited_memory=True
    )
    xtb_opt.optimize(mol=complex)
    complex.write(f'{name}_opt.mol')
    complex.dump(f'{name}_opt.json')

    return complex


def build_metal_organics(ligs):

    m4_FFs = molecule_building.metal_FFs(CN=4)
    m6_FFs = molecule_building.metal_FFs(CN=6)

    # Iterate over required library.
    for name in ligs:
        if ligs[name]['is_organic']:
            continue
        comp = ligs[name]
        opt_name = f'{name}_opt.mol'
        optjson_name = f'{name}_opt.json'
        if exists(optjson_name):
            continue
        print(f'.......building {name}')

        # Build metal atom.
        metal = molecule_building.build_metal(
            metal_smiles=comp['metal_smiles'],
            no_fgs=comp['no_metal_fgs']
        )
        # Get functional groups of binding atom to define coordination
        # sites.
        binding_atom = binding_atoms(string=comp['binding_atom'])
        binding_fgs = list(set([
            i.fg_type.name for i in binding_atom.func_groups
        ]))
        mc_topo = molecule_building.available_topologies(
            comp['metal_centre_topo']
        )
        metal_centre = molecule_building.build_metal_centre(
            metal=metal,
            topology=mc_topo,
            binding_atom=binding_atom,
            return_FG=binding_fgs
        )
        metal_centre.write(f'{name}_metal_centre.mol')

        # Load in organic BB.
        organic_BB = stk.BuildingBlock.init_from_file(
            f"{comp['organic_BB']}_opt.mol",
            functional_groups=comp['organic_FG']
        )
        # Build centre/complex using stk.
        ctopo = molecule_building.available_topologies(
            comp['ctopo']
        )
        n_metals = comp['no_metals']
        complex = stk.ConstructedMolecule(
            building_blocks=[metal_centre, organic_BB],
            topology_graph=ctopo,
            building_block_vertices={
                metal_centre: ctopo.vertices[:n_metals],
                organic_BB: ctopo.vertices[n_metals:]
            }
        )
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
        complex.dump(optjson_name)


def output_2d_images(ligs):

    if not exists('built_ligands/'):
        mkdir('built_ligands')

    # Draw 2D representation of all built molecules.
    mols = []
    names = []
    for name in ligs:
        comp = ligs[name]
        optjson_name = f'{name}_opt.json'
        if ligs[name]['is_organic']:
            BB = stk.BuildingBlock.load(optjson_name)
            MOL = rdkit.MolFromSmiles(
                rdkit.MolToSmiles(BB.to_rdkit_mol())
            )
            mols.append(MOL)
            AR = calculate_binding_AR(
                stk.BuildingBlock.init_from_molecule(
                    BB,
                    functional_groups=['bromine']
                )
            )
            label = f'{name}'
            if AR is not None:
                label = f'{label}: {round(AR,2)}'
            names.append(label)
        else:
            BB = stk.BuildingBlock.load(
                f"{comp['organic_BB']}_opt.json"
            )
            smi = rdkit.MolToSmiles(BB.to_rdkit_mol())
            msmi = comp['metal_smiles']
            tot_smi = f'{msmi}.{smi}'
            MOL = rdkit.MolFromSmiles(tot_smi)
            mols.append(MOL)
            top_str = comp['ctopo']
            label = f'{name}\n{top_str}'
            AR = calculate_binding_AR(
                stk.BuildingBlock.init_from_molecule(
                    BB,
                    functional_groups=['bromine']
                )
            )
            if AR is not None:
                label = f'{label}: {round(AR,2)}'
            names.append(label)

        atools.draw_mol_to_svg(
            mol=MOL,
            filename=f'built_ligands/{name}_opt.svg'
        )

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
