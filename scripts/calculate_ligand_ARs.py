#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot ligand AR as a function of defining atom.

Author: Andrew Tarzia

Date Created: 23 Feb 2021

"""

import sys
import numpy as np
from os.path import exists
from glob import glob
import matplotlib.pyplot as plt
import os
import shutil
from rdkit.Chem import AllChem as rdkit
import json

import stk

from utilities import (
    calculate_binding_AR,
    calculate_binding_ABS,
    get_planar_conformer,
    mol_list2grid,
    convert_lig_names_from_cage,
)
from cubeface import CubeFace
from facebuildingblock import FaceBuildingBlock


def run_reaction(ligand):

    corner_bb = stk.BuildingBlock(
        smiles='C1=CC=NC(=C1)C=NBr',
        functional_groups=[stk.BromoFactory()],
    )

    new_ligand = stk.ConstructedMolecule(
        topology_graph=CubeFace(
            building_blocks=(corner_bb, ligand),
        )
    )

    return new_ligand


def plot_parities(C_data, Br_data, N_data, output_name):

    fig, ax = plt.subplots(figsize=(8, 5))
    for lig in C_data:
        car = C_data[lig]
        brar = Br_data[lig]
        nar = N_data[lig]

        x = nar
        y1 = car
        y2 = brar
        ax.scatter(
            x,
            y1,
            c='gray',
            edgecolors='k',
            marker='o',
            alpha=1.0,
            s=120
        )
        ax.scatter(
            x,
            y2,
            c='sandybrown',
            edgecolors='k',
            marker='o',
            alpha=1.0,
            s=120
        )
        name = convert_lig_names_from_cage(lig)
        ax.text(x+0.03, y1, name, fontsize=12)
        ax.text(x+0.03, y2, name, fontsize=12)

    for c, name in zip(['gray', 'sandybrown'], ['C-C', 'Br-Br']):
        ax.scatter(
            -100, -100,
            c=c,
            edgecolors='k',
            marker='o',
            alpha=1.0,
            s=120,
            label=name,
        )


    ax.tick_params(axis='both', which='major', labelsize=16)
    if 'ARs' in output_name:
        ax.plot(
            np.linspace(1, 3), np.linspace(1, 3), '--', c='k', lw=2
        )
        ax.set_xlabel('N-N aspect ratio [1:X]', fontsize=16)
        ax.set_ylabel('test aspect ratio [1:X]', fontsize=16)
        ax.set_xlim(1.0, 2.6)
        ax.set_ylim(1.0, 2.6)
    elif 'ABs' in output_name:
        ax.plot(
            np.linspace(0, 8), np.linspace(0, 8), '--', c='k', lw=2
        )
        ax.set_xlabel('N-N difference [$\mathrm{\AA}$]', fontsize=16)
        ax.set_ylabel('test difference [$\mathrm{\AA}$]', fontsize=16)
        ax.set_xlim(0.0, 8)
        ax.set_ylim(0.0, 8)
    ax.legend(fontsize=16, ncol=2)
    fig.savefig(
        output_name,
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def write_values(C_data, Br_data, N_data, output_file):

    AR_file = output_file
    final_dict = {}
    for lig in C_data:
        final_dict[lig] = {
            'C': C_data[lig], 'N': N_data[lig], 'Br': Br_data[lig]
        }

    with open(AR_file, 'w') as f:
        json.dump(final_dict, f, indent=4)


def plot_values(C_data, Br_data, N_data, output_name):

    fig, ax = plt.subplots(figsize=(8, 5))
    cars = []
    nars = []
    brars = []
    names = []
    positions = []
    for i, lig in enumerate(C_data):
        cars.append(C_data[lig])
        brars.append(Br_data[lig])
        nars.append(N_data[lig])
        names.append(convert_lig_names_from_cage(lig))
        positions.append(convert_lig_names_from_cage(lig, as_int=True))

    # if 'ABs' not in output_name:
    #     ax.scatter(
    #         positions, cars,
    #         c='gray',
    #         edgecolors='k',
    #         marker='o',
    #         alpha=1.0,
    #         s=120,
    #         label='C-C',
    #     )

    #     ax.scatter(
    #         positions, nars,
    #         c='skyblue',
    #         edgecolors='k',
    #         marker='X',
    #         alpha=1.0,
    #         s=120,
    #         label='N-N',
    #     )

    ax.scatter(
        positions, brars,
        c='sandybrown',
        edgecolors='k',
        marker='s',
        alpha=1.0,
        s=120,
        label='Br-Br',
    )

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('subcomponent', fontsize=16)
    ax.set_xlim(0.0, max(positions)+1)
    if 'ARs' in output_name:
        ax.set_ylabel('aspect ratio [1:X]', fontsize=16)
        ax.set_ylim(1.0, 3)
    elif 'ABs' in output_name:
        ax.set_ylabel('difference [$\mathrm{\AA}$]', fontsize=16)
        ax.set_ylim(0.0, 8)
    ax.set_xticks(positions)
    ax.set_xticklabels(names)
    # ax.legend(fontsize=16, ncol=3)
    fig.savefig(
        output_name,
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def output_2d_images(C_data, Br_data, N_data):

    if os.path.exists('built_ligands'):
        shutil.rmtree('built_ligands')

    os.mkdir('built_ligands')

    # Draw 2D representation of all built molecules.
    mols = []
    cnames = []
    bnames = []
    nnames = []
    for name in C_data:
        car = C_data[name]
        brar = Br_data[name]
        nar = N_data[name]
        opt_name = f'{name}_opt.mol'
        mol = stk.BuildingBlock.init_from_file(opt_name).to_rdkit_mol()
        mol.RemoveAllConformers()
        mol = rdkit.RemoveHs(mol)
        rdkit.SanitizeMol(mol)
        mols.append(mol)
        label = f'{name}: {round(car, 2)}'
        cnames.append(label)
        label = f'{name}: {round(brar, 2)}'
        bnames.append(label)
        label = f'{name}: {round(nar, 2)}'
        nnames.append(label)

    # Draw 2D representation of all built molecules.
    mol_list2grid(
        molecules=mols,
        names=cnames,
        filename='built_ligands/cbuilt_ligands',
        mol_per_row=3,
        maxrows=3,
        subImgSize=(250, 200)
    )
    mol_list2grid(
        molecules=mols,
        names=bnames,
        filename='built_ligands/bbuilt_ligands',
        mol_per_row=3,
        maxrows=3,
        subImgSize=(250, 200)
    )
    mol_list2grid(
        molecules=mols,
        names=nnames,
        filename='built_ligands/nbuilt_ligands',
        mol_per_row=3,
        maxrows=3,
        subImgSize=(250, 200)
    )

    return


def main():
    first_line = (
        'Usage: calculate_ligand_ARs.py'
    )
    if (not len(sys.argv) == 1):
        print(f"""
{first_line}

    """)
        sys.exit()

    C_ars = {}
    Br_ars = {}
    N_ars = {}
    C_abs = {}
    Br_abs = {}
    N_abs = {}
    for l_file in glob('*_planar.mol'):
        ligand = l_file.replace('_planar.mol', '')
        planar_mol = FaceBuildingBlock.init_from_file(
            l_file, [stk.BromoFactory()],
        )
        if planar_mol.get_num_functional_groups() != 4:
            continue

        print(f'doing {ligand}')

        reacted_planar_file = f'{ligand}_reacted_planar.mol'
        if exists(reacted_planar_file):
            reacted_planar_mol = stk.BuildingBlock.init_from_file(
                reacted_planar_file
            )
        else:
            print(f'getting planar, reacted {ligand}')
            reacted_mol = run_reaction(planar_mol)
            reacted_planar_mol = get_planar_conformer(reacted_mol, 300)
            reacted_planar_mol.write(reacted_planar_file)

        reacted_planar_mol = stk.BuildingBlock.init_from_molecule(
            molecule=reacted_planar_mol,
            functional_groups=(
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]-[#7X2]=[#6X3H1]-[#6X3!H1]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ),
        )
        if reacted_planar_mol.get_num_functional_groups() != 4:
            raise ValueError('ligand reacted did not have 4 FGs')

        car = calculate_binding_AR(planar_mol, atom_ids=None)
        cab = calculate_binding_ABS(planar_mol, atom_ids=None)
        br_atom_ids = [
            list(fg.get_deleter_ids())
            for fg in planar_mol.get_functional_groups()
        ]
        brar = calculate_binding_AR(planar_mol, atom_ids=br_atom_ids)
        brab = calculate_binding_ABS(planar_mol, atom_ids=br_atom_ids)
        n_atom_ids = [
            list(fg.get_bonder_ids())
            for fg in reacted_planar_mol.get_functional_groups()
        ]
        nar = calculate_binding_AR(
            mol=reacted_planar_mol, atom_ids=n_atom_ids
        )
        nab = calculate_binding_ABS(
            mol=reacted_planar_mol, atom_ids=n_atom_ids
        )
        C_ars[ligand] = car
        C_abs[ligand] = cab
        Br_ars[ligand] = brar
        Br_abs[ligand] = brab
        N_ars[ligand] = nar
        N_abs[ligand] = nab

    plot_parities(C_abs, Br_abs, N_abs, 'all_ABs.pdf')
    plot_values(C_abs, Br_abs, N_abs, 'all_ABs_bar.pdf')
    write_values(C_abs, Br_abs, N_abs, 'ligand_ABs.json')

    plot_parities(C_ars, Br_ars, N_ars, 'all_ARs.pdf')
    plot_values(C_ars, Br_ars, N_ars, 'all_ARs_bar.pdf')
    write_values(C_ars, Br_ars, N_ars, 'ligand_ARs.json')

    output_2d_images(C_abs, Br_abs, N_abs)


if __name__ == "__main__":
    main()
