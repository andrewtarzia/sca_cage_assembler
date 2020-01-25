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
        output = f'{name}_opt.mol'
        if exists(output):
            continue
        print(f'doing {name}')
        mol = stk.BuildingBlock(smiles=smi)
        optimizer.optimize(mol)
        mol.write(output)

    return


def metal_containing_ligands():
    """
    Defines how to build metal containing ligands.

    Uses smiles strings defined in ligand lib and stk topologies.

    """

    m_FFs =  {
        30: 'Zn4+2', 28: 'Ni4+2',
        78: 'Pt4+2', 46: 'Pd4+2',
        45: 'Rh6+3', 42: 'Mo4f2'
    }

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


def build_metal_organics(metal_lig_lib, ligs):

    return


def output_2d_image():
    # Draw 2D representation of all built molecules.
    mol_list = []
    name_list = []
    opt_mols = sorted(glob.glob('_opt.mol'))
    for mol in opt_mols:
        name_list.append(mol.replace('_opt.mol'))
        MOL = rdkit.MolFromMolFile(mol)
        mol_list.append(MOL)

    atools.mol_list2grid(
        molecules=mol_list,
        filename='built_ligands',
        names=name_list,
        mol_per_row=3,
        maxrows=3,
        subImgSize=(200, 200)
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
    output_2d_image()
    sys.exit()


if __name__ == "__main__":
    main()
