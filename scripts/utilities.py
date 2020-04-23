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
from rdkit.Chem import AllChem as rdkit

import atools


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


def build_conformers(mol, N):
    molecule = mol.to_rdkit_mol()
    molecule.RemoveAllConformers()

    cids = rdkit.EmbedMultipleConfs(
        mol=molecule,
        numConfs=N,
        randomSeed=1000,
        useExpTorsionAnglePrefs=True,
        useBasicKnowledge=True
    )
    print(f'there are {molecule.GetNumConformers()} conformers')
    input()
    return cids, molecule


def get_lowest_energy_conformer(
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
    Get lowest energy conformer of molecule.

    Method:
        1) ETKDG conformer search on molecule
        2) xTB normal optimisation of each conformer
        3) xTB opt_level optimisation of lowest energy conformer
        4) save file

    """

    # Run ETKDG on molecule.
    print(f'....running ETKDG on {name}')
    cids, confs = build_conformers(mol, N=100)

    # Optimize all conformers at normal level with xTB.
    low_e_conf_id = -100
    low_e = 10E20
    for cid in cids:
        name_ = f'{name}_confs/c_{cid}'
        ey_file = f'{name}_confs/c_{cid}_eyout'
        print(name, ey_file)
        mol = atools.update_from_rdkit_conf(
            mol,
            confs,
            conf_id=cid
        )
        mol.write(f'temp_c_{cid}.mol')

        # Optimize.
        opt_mol = optimize_molecule(
            name=name_+'_opt',
            mol=mol,
            opt_level='normal'
        )

        # Get energy.
        calculate_energy(
            name=name_+'_ey',
            mol=opt_mol,
            ey_file=ey_file
        )
        ey = read_ey(ey_file)
        if ey < low_e:
            low_e_conf_id = cid
            low_e = ey
        print(ey, low_e, low_e_conf_id, cid)

    # Get lowest energy conformer.
    low_e_conf = atools.update_from_rdkit_conf(
        mol,
        confs,
        conf_id=low_e_conf_id
    )
    low_e_conf.write('temp_pre_opt.mol')

    # Optimize lowest energy conformer at opt_level.
    low_e_conf = optimize_molecule(
        name=name_+'low_e_opt',
        mol=low_e_conf,
        opt_level=opt_level
    )
    low_e_conf.write('temp_post_opt.mol')
    print(low_e_conf_id)

    # Return molecule.
    return low_e_conf


def calculate_binding_AR(mol):
    """
    Calculate ligand aspect ratio based on binder positions.

    Defined as:

    """

    print(mol.func_groups)
    fg_ids = range(len(mol.func_groups))
    print(fg_ids)
    if len(mol.func_groups) != 4:
        return None

    binder_atom_ids = [
        list(mol.get_bonder_ids(fg_ids=[i]))
        for i in fg_ids
    ]
    print(binder_atom_ids)
    binder_atom_dists = sorted(
        list(mol.get_bonder_distances()),
        key=lambda a: a[2]
    )

    print(binder_atom_dists)
    far_binder_pair = (
        binder_atom_dists[-1][0],
        binder_atom_dists[-1][1]
    )
    print(far_binder_pair)
    ARs = []
    for fg_id in far_binder_pair:
        print(fg_id)
        ds = sorted([
            i[2] for i in binder_atom_dists
            if fg_id in (i[0], i[1])
        ])
        print(ds)
        AR = ds[1]/min(ds)
        print(AR)
        ARs.append(AR)

    ligand_AR = sum(ARs)/len(ARs)
    print(ARs, ligand_AR)
    return ligand_AR
