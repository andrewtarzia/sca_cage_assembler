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
import matplotlib.pyplot as plt
import json

import stk

from utilities import calculate_binding_AR, get_planar_conformer
from cubeface import CubeFace


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


def plot_AR_parities(C_ars, Br_ars, N_ars):

    fig, ax = plt.subplots(figsize=(8, 5))
    for lig in C_ars:
        car = C_ars[lig]
        brar = Br_ars[lig]
        nar = N_ars[lig]

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
        name = lig.replace('quad2_', '')
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

    ax.plot(np.linspace(1, 3), np.linspace(1, 3), '--', c='k', lw=2)

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('N-N aspect ratio [1:X]', fontsize=16)
    ax.set_ylabel('test aspect ratio [1:X]', fontsize=16)
    ax.set_xlim(1.0, 2.6)
    ax.set_ylim(1.0, 2.6)
    ax.legend(fontsize=16, ncol=2)
    fig.savefig(
        'all_ARs.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def write_ARs(C_ars, Br_ars, N_ars):

    AR_file = 'ligand_ARs.json'
    final_dict = {}
    for lig in C_ars:
        final_dict[lig] = {
            'C': C_ars[lig], 'N': N_ars[lig], 'Br': Br_ars[lig]
        }

    with open(AR_file, 'w') as f:
        json.dump(final_dict, f, indent=4)


def plot_ARs(C_ars, Br_ars, N_ars):

    fig, ax = plt.subplots(figsize=(8, 5))
    cars = []
    nars = []
    brars = []
    names = []
    positions = []
    for i, lig in enumerate(C_ars):
        cars.append(C_ars[lig])
        brars.append(Br_ars[lig])
        nars.append(N_ars[lig])
        name = lig.replace('quad2_', '')
        names.append(name)
        positions.append(i+1)

    ax.scatter(
        positions, cars,
        c='gray',
        edgecolors='k',
        marker='o',
        alpha=1.0,
        s=120,
        label='C-C',
    )

    ax.scatter(
        positions, nars,
        c='skyblue',
        edgecolors='k',
        marker='X',
        alpha=1.0,
        s=120,
        label='N-N',
    )

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
    ax.set_xlabel('ligand name', fontsize=16)
    ax.set_ylabel('aspect ratio [1:X]', fontsize=16)
    ax.set_xlim(0.0, max(positions)+1)
    ax.set_ylim(1.0, 3)
    ax.set_xticks(positions)
    ax.set_xticklabels(names)
    ax.legend(fontsize=16, ncol=3)
    fig.savefig(
        'all_ARs_bar.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def main():
    first_line = (
        'Usage: plot_AR_comparisons.py'
    )
    if (not len(sys.argv) == 1):
        print(f"""
{first_line}

    """)
        sys.exit()

    # Skip ligands not in database of 10.
    dataset_of_10 = [
        'quad2_1', 'quad2_12', 'quad2_2', 'quad2_3', 'quad2_8',
        'quad2_9', 'quad2_10', 'quad2_5', 'quad2_16', 'quad2_17',
    ]

    C_ars = {}
    Br_ars = {}
    N_ars = {}
    for ligand in dataset_of_10:
        print(ligand)
        planar_mol = stk.BuildingBlock.init_from_file(
            f'{ligand}_planar.mol', [stk.BromoFactory()],
        )
        print(planar_mol)

        reacted_planar_file = f'{ligand}_reacted_planar.mol'
        if exists(reacted_planar_file):
            reacted_planar_mol = stk.BuildingBlock.init_from_file(
                reacted_planar_file
            )
        else:
            print(f'running {ligand}')
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

        print(reacted_planar_mol)

        car = calculate_binding_AR(planar_mol, atom_ids=None)
        br_atom_ids = [
            list(fg.get_deleter_ids())
            for fg in planar_mol.get_functional_groups()
        ]
        brar = calculate_binding_AR(planar_mol, atom_ids=br_atom_ids)
        n_atom_ids = [
            list(fg.get_bonder_ids())
            for fg in reacted_planar_mol.get_functional_groups()
        ]
        nar = calculate_binding_AR(
            mol=reacted_planar_mol, atom_ids=n_atom_ids
        )
        print(car)
        print(brar)
        print(nar)
        C_ars[ligand] = car
        Br_ars[ligand] = brar
        N_ars[ligand] = nar

    plot_AR_parities(C_ars, Br_ars, N_ars)
    plot_ARs(C_ars, Br_ars, N_ars)
    write_ARs(C_ars, Br_ars, N_ars)


if __name__ == "__main__":
    main()
