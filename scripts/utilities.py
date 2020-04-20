#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Utilities module.

Author: Andrew Tarzia

Date Created: 15 Mar 2020
"""

import json
import stk


def read_lib(lib_file):
    """
    Read lib file.

    Returns dictionary.

    """

    print(f'reading {lib_file}')
    with open(lib_file, 'rb') as f:
        lib = json.load(f)

    return lib


def read_ey(file):
    """
    Read the energy (kJ/mol from GFN2-xTB) from a .ey file.

    """

    with open(file, 'r') as f:
        lines = f.readlines()
        ey = float(lines[0].rstrip())

    return ey*2625.5


def calculate_energy(
    name,
    mol,
    ey_file,
    charge=0,
    no_unpaired_e=0,
    solvent=None
):
    """
    Calculate GFN-xTB energy of molecule.

    """
    print(f'....getting energy of {name}')
    if solvent is None:
        solvent_str = None
        solvent_grid = 'normal'
    else:
        solvent_str, solvent_grid = solvent
    xtb_energy = stk.XTBEnergy(
        xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
        output_dir=f'{name}_ey',
        num_cores=6,
        charge=charge,
        num_unpaired_electrons=no_unpaired_e,
        electronic_temperature=300,
        unlimited_memory=True,
        calculate_free_energy=False,
        solvent=solvent_str,
        solvent_grid=solvent_grid
    )
    energy = xtb_energy.get_energy(mol)

    with open(ey_file, 'w') as f:
        f.write(str(energy))


def optimize_molecule(
    name,
    mol,
    opt_level='extreme',
    charge=0,
    no_unpaired_e=0,
    max_runs=1,
    calc_hessian=False,
    solvent=None
):
    """
    Run simple GFN-xTB optimisation of molecule.

    """
    print(f'....optimizing {name}')
    if solvent is None:
        solvent_str = None
        solvent_grid = 'normal'
    else:
        solvent_str, solvent_grid = solvent
    xtb_opt = stk.XTB(
        xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
        output_dir=f'{name}_opt',
        gfn_version=2,
        num_cores=6,
        opt_level=opt_level,
        charge=charge,
        num_unpaired_electrons=no_unpaired_e,
        max_runs=max_runs,
        calculate_hessian=calc_hessian,
        unlimited_memory=True,
        solvent=solvent_str,
        solvent_grid=solvent_grid
    )
    xtb_opt.optimize(mol=mol)

    return mol
