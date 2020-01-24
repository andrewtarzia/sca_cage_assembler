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

    # output_name : (
    #     organic ligand name,
    #     metal centre topology_graph,
    #     metal smiles,
    #     no. FGs
    # )

    m_ligands = {
        'quad4_3': (
            'quad4_prec_4',
            stk.metal_centre.Porphyrin,
            '[Ni+2]',
            4
        ),
        'quad4_4': (
            'quad4_prec_5',
            stk.metal_centre.Porphyrin,
            '[Ni+2]',
            4
        ),
        'quad4_5': (
            'quad4_prec_4',
            stk.metal_centre.Porphyrin,
            '[Zn+2]',
            4
        ),
        'quad4_6': (
            'quad4_prec_5',
            stk.metal_centre.Porphyrin,
            '[Zn+2]',
            4
        ),
        'quad4_7': (
            'quad4_prec_1',
            stk.metal_centre.SquarePlanar,
            '[Pt+2]',
            4
        ),
        'quad4_10': (
            'quad4_prec_1',
            stk.metal_centre.SquarePlanar,
            '[Pd+2]',
            4
        ),
        'quad4_8': (
            'quad4_prec_2',
            stk.metal_centre.Paddlewheel,
            '[Rh+2]',
            4
        ),
        'quad4_9': (
            'quad4_prec_2',
            stk.metal_centre.Paddlewheel,
            '[Mo+2]',
            4
        ),
        'quad4_11': (
            'quad4_prec_3',
            stk.metal_centre.SquarePlanar,
            '[Pt+2]',
            4
        ),
        'quad4_12': (
            'quad4_prec_3',
            stk.metal_centre.SquarePlanar,
            '[Pd+2]',
            4
        ),
    }

    return m_ligands


def build_metal_organics(metal_lig_lib, ligs):

    return


def output_2d_image():

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
