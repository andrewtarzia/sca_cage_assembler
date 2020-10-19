#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build and analyse faces from ligand library.

Author: Andrew Tarzia

Date Created: 19 Oct 2020

"""

import sys
from os.path import exists, join
from glob import glob
import numpy as np

import stk
import stko

from molecule_building import metal_FFs
def load_complex(filename):

    fgfactory = stk.SmartsFunctionalGroupFactory(
        smarts='[#7X3]~[#6]~[#6]~[#7X3]~[#35]',
        bonders=(3, ),
        deleters=(4, ),
        placers=(0, 1, 2, 3),
    )

    name = filename.replace('.mol', '')
    # Need to define more than one placer id for complexes to ensure
    # alignment -- use the NCCN plane (attached to the Br) to define
    # the orientation of the complex.
    complex = stk.BuildingBlock.init_from_file(
        filename,
        functional_groups=[fgfactory]
    )

    return name, complex


def optimize_complex(complex, name):

    opt_name = f'{name}_opt.mol'
    if exists(opt_name):
        return complex.with_structure_from_file(opt_name)
    else:
        print(f'doing UFF4MOF optimisation for {name}')
        gulp_opt = stko.GulpUFFOptimizer(
            gulp_path='/home/atarzia/software/gulp-5.1/Src/gulp/gulp',
            metal_FF=metal_FFs(CN=6),
            output_dir=f'{name}_uff1'
        )
        gulp_opt.assign_FF(complex)
        complex = gulp_opt.optimize(mol=complex)
        complex.write(f'{name}_uff1.mol')

        print(f'doing xTB optimisation for {name}')
        xtb_opt = stko.XTB(
            xtb_path='/home/atarzia/software/xtb-6.3.1/bin/xtb',
            output_dir=f'{name}_xtb',
            gfn_version=2,
            num_cores=6,
            opt_level='tight',
            charge=2,
            num_unpaired_electrons=0,
            max_runs=1,
            calculate_hessian=False,
            unlimited_memory=True
        )
        complex = xtb_opt.optimize(mol=complex)
        complex.write(f'{name}_opt.mol')


def load_ligands(directory):

    ligands = {}
    for lig in glob(join(directory, f'quad2*_opt.mol')):
        l_name = lig.replace(directory, '').replace('_opt.mol', '')
        ligands[l_name] = stk.BuildingBlock.init_from_file(
            lig,
            functional_groups=[stk.BromoFactory()],
        )

    return ligands

def main():
    first_line = (
        'Usage: face_analysis.py '
        'lig_directory'
    )
    if (not len(sys.argv) == 2):
        print(f"""
{first_line}

    ligand_directory : (str)
        Directory with required ligand structures.

    """)
        sys.exit()
    else:
        ligand_directory = sys.argv[1]

    # Define two metal building blocks (lambda, delta).
    del_name, del_complex = load_complex('cl1_zn_oct_del_face.mol')
    lam_name, lam_complex = load_complex('cl1_zn_oct_lam_face.mol')

    # Optimise both complexes.
    del_complex = optimize_complex(del_complex, del_name)
    lam_complex = optimize_complex(lam_complex, lam_name)

    # Load in each ligand structure.
    ligands = load_ligands(ligand_directory)


if __name__ == '__main__':
    main()
