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
import glob
from rdkit.Chem import AllChem as rdkit
import stk

import atools
import Building


def read_lig_lib(lib_file):
    """
    Read ligand lib file.

    Returns dictionary of format:

    ligs[name] = (smiles, flag)

    """
    ligs = {}
    with open(lib_file, 'r') as f:
        lines = f.readlines()

    for line in lines[1:]:
        row = line.rstrip().split(',')
        ligs[row[0]] = (row[1], row[2])
    return ligs


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
        smi = ligs[name][0]
        input = f'manual/{name}.mol'
        output = f'{name}_opt.mol'
        jsonoutput = f'{name}_opt.json'
        if exists(jsonoutput):
            continue
        print(f'doing {name}')
        if exists(input):
            mol = stk.BuildingBlock.init_from_file(input)
        else:
            mol = stk.BuildingBlock(smiles=smi)
        optimizer.optimize(mol)
        mol.write(output)
        mol.dump(jsonoutput)

    return


def metal_containing_ligands():
    """
    Defines how to build metal containing ligands.

    Uses smiles strings defined in ligand lib and stk topologies.

    """

    m_FFs = Building.metal_FFs()

    m_ligands = {
        'quad4_3': {
            'organic_BB': 'quad4_prec_4',
            'organic_FG': ['CNC_metal'],
            'ctopo': stk.cage.Porphyrin(),
            'metal_centre_topo': stk.metal_centre.SquarePlanar(),
            'metal_smiles': '[Ni+2]',
            'no_metal_fgs': 4,
            'no_metals': 1,
            'net_charge': 0,
            'binding_atom': Building.build_atom(
                'N',
                FG='metal_bound_N'
            ),
            'metal_FF': m_FFs
        },
        'quad4_4': {
            'organic_BB': 'quad4_prec_5',
            'organic_FG': ['CNC_metal'],
            'ctopo': stk.cage.Porphyrin(),
            'metal_centre_topo': stk.metal_centre.SquarePlanar(),
            'metal_smiles': '[Ni+2]',
            'no_metal_fgs': 4,
            'no_metals': 1,
            'net_charge': 0,
            'binding_atom': Building.build_atom(
                'N',
                FG='metal_bound_N'
            ),
            'metal_FF': m_FFs
        },
        'quad4_5': {
            'organic_BB': 'quad4_prec_4',
            'organic_FG': ['CNC_metal'],
            'ctopo': stk.cage.Porphyrin(),
            'metal_centre_topo': stk.metal_centre.SquarePlanar(),
            'metal_smiles': '[Zn+2]',
            'no_metal_fgs': 4,
            'no_metals': 1,
            'net_charge': 0,
            'binding_atom': Building.build_atom(
                'N',
                FG='metal_bound_N'
            ),
            'metal_FF': m_FFs
        },
        'quad4_6': {
            'organic_BB': 'quad4_prec_5',
            'organic_FG': ['CNC_metal'],
            'ctopo': stk.cage.Porphyrin(),
            'metal_centre_topo': stk.metal_centre.SquarePlanar(),
            'metal_smiles': '[Zn+2]',
            'no_metal_fgs': 4,
            'no_metals': 1,
            'net_charge': 0,
            'binding_atom': Building.build_atom(
                'N',
                FG='metal_bound_N'
            ),
            'metal_FF': m_FFs
        },
        'quad4_7': {
            'organic_BB': 'quad4_prec_1',
            'organic_FG': ['pyridine_N_metal'],
            'ctopo': stk.cage.SquarePlanarMonodentate(),
            'metal_centre_topo': stk.metal_centre.SquarePlanar(),
            'metal_smiles': '[Pt+2]',
            'no_metal_fgs': 4,
            'no_metals': 1,
            'net_charge': 2,
            'binding_atom': Building.build_atom(
                'N',
                FG='metal_bound_N'
            ),
            'metal_FF': m_FFs
        },
        'quad4_10': {
            'organic_BB': 'quad4_prec_1',
            'organic_FG': ['pyridine_N_metal'],
            'ctopo': stk.cage.SquarePlanarMonodentate(),
            'metal_centre_topo': stk.metal_centre.SquarePlanar(),
            'metal_smiles': '[Pd+2]',
            'no_metal_fgs': 4,
            'no_metals': 1,
            'net_charge': 2,
            'binding_atom': Building.build_atom(
                'N',
                FG='metal_bound_N'
            ),
            'metal_FF': m_FFs
        },
        'quad4_8': {
            'organic_BB': 'quad4_prec_2',
            'organic_FG': ['CO_metal', 'COH_metal'],
            'ctopo': stk.cage.Paddlewheel(),
            'metal_centre_topo': stk.metal_centre.SquarePlanar(),
            'metal_smiles': '[Rh+2]',
            'no_metal_fgs': 4,
            'no_metals': 2,
            'net_charge': 0,
            'binding_atom': Building.build_atom(
                'O',
                FG='metal_bound_O'
            ),
            'metal_FF': m_FFs
        },
        'quad4_9': {
            'organic_BB': 'quad4_prec_2',
            'organic_FG': ['CO_metal', 'COH_metal'],
            'ctopo': stk.cage.Paddlewheel(),
            'metal_centre_topo': stk.metal_centre.SquarePlanar(),
            'metal_smiles': '[Mo+2]',
            'no_metal_fgs': 4,
            'no_metals': 2,
            'net_charge': 0,
            'binding_atom': Building.build_atom(
                'O',
                FG='metal_bound_O'
            ),
            'metal_FF': m_FFs
        },
        'quad4_11': {
            'organic_BB': 'quad4_prec_3',
            'organic_FG': ['pyridine_N_metal'],
            'ctopo': stk.cage.SquarePlanarMonodentate(),
            'metal_centre_topo': stk.metal_centre.SquarePlanar(),
            'metal_smiles': '[Pt+2]',
            'no_metal_fgs': 4,
            'no_metals': 1,
            'net_charge': 2,
            'binding_atom': Building.build_atom(
                'N',
                FG='metal_bound_N'
            ),
            'metal_FF': m_FFs
        },
        'quad4_12': {
            'organic_BB': 'quad4_prec_3',
            'organic_FG': ['pyridine_N_metal'],
            'ctopo': stk.cage.SquarePlanarMonodentate(),
            'metal_centre_topo': stk.metal_centre.SquarePlanar(),
            'metal_smiles': '[Pd+2]',
            'no_metal_fgs': 4,
            'no_metals': 1,
            'net_charge': 2,
            'binding_atom': Building.build_atom(
                'N',
                FG='metal_bound_N'
            ),
            'metal_FF': m_FFs
        },
    }

    return m_ligands


def optimise_metal_centre(name, charge, complex, metal_FF):

    print(f'doing UFF4MOF optimisation of {name}')
    gulp_opt = stk.GulpMetalOptimizer(
        gulp_path='/home/atarzia/software/gulp-5.1/Src/gulp/gulp',
        metal_FF=metal_FF,
        output_dir=f'{name}_uff1'
    )
    gulp_opt.assign_FF(complex)
    gulp_opt.optimize(mol=complex)
    complex.write(f'{name}_uff1.mol')
    complex.dump(f'{name}_uff1.json')

    print(f'doing XTB optimisation of {name}')
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


def build_metal_organics(metal_lig_lib, ligs):

    # Iterate over required metal-organic library.
    for name in metal_lig_lib:
        opt_name = f'{name}_opt.mol'
        optjson_name = f'{name}_opt.json'
        if exists(optjson_name):
            continue
        print(f'building {name}')
        comp = metal_lig_lib[name]

        # Build metal atom.
        metal = Building.build_metal(
            metal_smiles=comp['metal_smiles'],
            no_fgs=comp['no_metal_fgs']
        )
        # Get functional groups of binding atom to define coordination
        # sites.
        binding_fgs = list(set([
            i.fg_type.name for i in comp['binding_atom'].func_groups
        ]))
        metal_centre = Building.build_metal_centre(
            metal=metal,
            topology=comp['metal_centre_topo'],
            binding_atom=comp['binding_atom'],
            return_FG=binding_fgs
        )
        metal_centre.write(f'{name}_metal_centre.mol')

        # Load in organic BB.
        organic_BB = stk.BuildingBlock.init_from_file(
            f"{comp['organic_BB']}_opt.mol",
            functional_groups=comp['organic_FG']
        )
        # Build centre/complex using stk.
        ctopo = comp['ctopo']
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
        complex = optimise_metal_centre(
            name=name,
            charge=comp['net_charge'],
            complex=complex,
            metal_FF=comp['metal_FF']
        )
        # Output.
        complex.write(opt_name)
        complex.dump(optjson_name)

    return


def output_2d_images():
    # Draw 2D representation of all built molecules.
    opt_mols = sorted(glob.glob('*_opt.json'))
    for mol in opt_mols:
        name = mol.replace('_opt.json', '')
        MOL = stk.BuildingBlock.load(mol).to_rdkit_mol()
        rdkit.Compute2DCoords(MOL)
        atools.draw_mol_to_svg(
            mol=MOL,
            filename=f'built_ligands/{name}_opt.svg'
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

    print(f'reading {lib_file}')
    ligs = read_lig_lib(lib_file)

    # Build and optimise all organic molecules in lib.
    build_organics(ligs)

    # Build and optimise all metal containing ligands.
    metal_lig_lib = metal_containing_ligands()
    build_metal_organics(metal_lig_lib, ligs)

    # Produce image of all built molecules.
    output_2d_images()
    sys.exit()


if __name__ == "__main__":
    main()
