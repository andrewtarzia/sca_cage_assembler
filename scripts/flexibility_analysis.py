#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to calculate flexibility measure of ligands.

Author: Andrew Tarzia

Date Created: 21 Oct 2020

"""

import numpy as np
from itertools import combinations
import sys
import os
import matplotlib.pyplot as plt
import json

import stk

from plotting import colors_i_like, histogram_plot_N
from utilities import (
    split_xyz_file,
    get_atom_distance,
    calculate_molecule_planarity,
    get_lowest_energy_conformer,
    read_lib,
)
import env_set


def get_crest_ensemble_data(crest_directory):

    with open(f'{crest_directory}/crest.output', 'r') as f:
        for line in f.readlines():
            # Get number of conformers.
            if ' number of unique conformers for further calc' in line:
                no_conformers = (
                    int(line.rstrip().split(' ')[-1])
                )

            # Get number of rotamers.
            if 'total number unique points considered further' in line:
                no_rotamers = (
                    int(line.rstrip().split(' ')[-1])
                )

    return no_rotamers, no_conformers


def is_single_binder(molecule):
    for fg in molecule.get_functional_groups():
        if len(list(fg.get_bonder_ids())) > 1:
            return False
    return True

def calculate_long_axis_distance(molecule, conformer_files):

    if not is_single_binder(molecule):
        raise ValueError(
            f'{molecule} has FGs with more than one binder.'
        )
    long_axis_atom_pair = get_long_axis_atoms(molecule)
    pair1_cents = [
        molecule.with_structure_from_file(i).get_centroid(
            atom_ids=tuple(long_axis_atom_pair[0])
        )
        for i in conformer_files
    ]
    pair2_cents = [
        molecule.with_structure_from_file(i).get_centroid(
            atom_ids=tuple(long_axis_atom_pair[1])
        )
        for i in conformer_files
    ]

    la_dist = [
        np.linalg.norm(i-j)
        for i, j in zip(pair1_cents, pair2_cents)
    ]
    return la_dist


def plot_long_axis_deviation(measures, name, crest=False):
    # Can assume the first one in the list of measures is the lowest
    # energy conformer.
    fig, ax = histogram_plot_N(
        Y=[i-measures[0] for i in measures],
        X_range=(-1.5, 1.5),
        width=0.05,
        alpha=1.0,
        color=colors_i_like()[1],
        edgecolor=colors_i_like()[1],
        xtitle=r'long axis deviation [$\mathrm{\AA}$]',
        N=1
    )
    range = abs(max(measures) - min(measures))
    ax.set_title(f'range = {round(range, 2)}', fontsize=16)
    fig.tight_layout()
    if crest:
        filename = f'{name}_lapC_dist.pdf'
    else:
        filename = f'{name}_lap_dist.pdf'
    fig.savefig(
        filename,
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_plane_deviation(measures, name, crest=False):

    fig, ax = histogram_plot_N(
        Y=measures,
        X_range=(0, 200),
        width=4,
        alpha=1.0,
        color=colors_i_like()[1],
        edgecolor=colors_i_like()[1],
        xtitle=r'plane deviation [$\mathrm{\AA}$]',
        N=1
    )
    range = abs(max(measures) - min(measures))
    ax.set_title(f'range = {round(range, 2)}', fontsize=16)
    fig.tight_layout()
    if crest:
        filename = f'{name}_AAplanedevC_dist.pdf'
    else:
        filename = f'{name}_AAplanedev_dist.pdf'
    fig.savefig(
        filename,
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_binder_plane_deviation(measures, name, crest=False):

    fig, ax = histogram_plot_N(
        Y=measures,
        X_range=(0, 30),
        width=0.2,
        alpha=1.0,
        color=colors_i_like()[1],
        edgecolor=colors_i_like()[1],
        xtitle=r'binder plane deviation [$\mathrm{\AA}$]',
        N=1
    )
    range = abs(max(measures) - min(measures))
    ax.set_title(f'range = {round(range, 2)}', fontsize=16)
    fig.tight_layout()
    if crest:
        filename = f'{name}_planedevC_dist.pdf'
    else:
        filename = f'{name}_planedev_dist.pdf'
    fig.savefig(
        filename,
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def get_long_axis_atoms(molecule):

    binder_atom_ids = [
        list(fg.get_bonder_ids())
        for fg in molecule.get_functional_groups()
    ]
    binder_atom_dists = sorted(
        [
            (idx1, idx2, get_atom_distance(
                molecule,
                idx1,
                idx2
            ))
            for idx1, idx2 in combinations(binder_atom_ids, r=2)
        ],
        key=lambda a: a[2]
    )
    # Can assume the ordering of the binder atom distances:
    # 0, 1: short vectors
    # 2, 3: long vectors
    # 4, 5: diagonal vectors
    # This fails when the molecule is not sufficiently anisotropic,
    # at which point it will not matter.
    short_vector_fg_1 = (
        binder_atom_dists[0][0], binder_atom_dists[0][1]
    )
    short_vector_fg_2 = (
        (binder_atom_dists[1][0], binder_atom_dists[1][1])
        if (
            binder_atom_dists[1][0] not in short_vector_fg_1 and
            binder_atom_dists[1][1] not in short_vector_fg_1
        ) else
        (binder_atom_dists[2][0], binder_atom_dists[2][1])
    )

    long_axis_atom_pairs = (
        [i[0] for i in short_vector_fg_1],
        [i[0] for i in short_vector_fg_2]
    )

    return long_axis_atom_pairs


def main():
    first_line = (
        'Usage: flexibility_analysis.py ligand_directory '
        'ligand_lib_file'
    )
    if (not len(sys.argv) == 3):
        print(f"""
{first_line}

    ligand_directory : (str)
        Directory with required ligand structures.

    ligand_lib_file : (str)
        File containing ligand information (XXXXX)

    """)
        sys.exit()
    else:
        ligand_directory = sys.argv[1]
        ligand_lib_file = sys.argv[2]

    # Load in each ligand structure.
    ligand_lib = read_lib(ligand_lib_file)
    ligand_structures = {}

    for name in ligand_lib:
        structure_file = os.path.join(
            ligand_directory, f'{name}_opt.mol'
        )
        bb = stk.BuildingBlock.init_from_file(
            structure_file,
            functional_groups=[stk.BromoFactory()],
        )
        if ligand_lib[name]['calculate_flex']:
            ligand_structures[name] = bb

    print(f'there are {len(ligand_structures)} structures\n')
    for name in ligand_structures:
        lig_structure = ligand_structures[name]
        crest_output_file = f'{name}_flex_measure.json'
        crest_data = {}

        low_e_conformer_output = f'../{name}_loweconf.mol'
        conf_dir = f'{name}_xtbcrest_confs'
        if not os.path.exists(conf_dir):
            os.mkdir(conf_dir)
        # Crest part.
        if not os.path.exists(low_e_conformer_output):
            new_molecule = get_lowest_energy_conformer(
                name=name,
                mol=lig_structure,
                conf_dir=conf_dir,
                settings=env_set.crest_conformer_settings(
                    solvent=None,
                ),
            )
            new_molecule.write(low_e_conformer_output)

        # Extract some measure of conformer ensemble size.
        no_rotamers, no_conformers = get_crest_ensemble_data(
            crest_directory=f'{conf_dir}'
        )
        crest_data['no_rotamers'] = no_rotamers
        crest_data['no_conformers'] = no_conformers

        # Analyse all conformers from CREST.
        crest_conformer_files = split_xyz_file(
            num_atoms=lig_structure.get_num_atoms(),
            xyz_file=f'{conf_dir}/crest_conformers.xyz',
        )
        print(f'{name} has {len(crest_conformer_files)} conformers')
        # Plane deviations.
        crest_data['plane_deviations'] = [
            calculate_molecule_planarity(
                mol=lig_structure.with_structure_from_file(i),
            )
            for i in crest_conformer_files
        ]

        print(
            f": {name}, num rotamers: {crest_data['no_rotamers']}, "
            f"num conformers: {crest_data['no_conformers']}"
        )
        plot_plane_deviation(
            measures=crest_data['plane_deviations'],
            name=name,
            crest=True,
        )

        # Do functional group dependant analysis.
        if lig_structure.get_num_functional_groups() == 4:
            # Tetratopic, single binder specific measures.
            long_axis_distances = calculate_long_axis_distance(
                molecule=lig_structure,
                conformer_files=crest_conformer_files,
            )
            dist_width = abs(
                max(long_axis_distances)
                -min(long_axis_distances)
            )
            print(f':: {name}, dist width = {dist_width}')
            plot_long_axis_deviation(
                measures=long_axis_distances,
                name=name,
                crest=True,
            )
            crest_data['long_axis_distances'] = long_axis_distances
        if lig_structure.get_num_functional_groups() == 3:
            # Tritopic specific measures.
            pass
        if lig_structure.get_num_functional_groups() == 2:
            # Ditopic specific measures.
            pass
        if lig_structure.get_num_functional_groups() > 1:
            # At least two functional groups.
            binder_ids = [
                fg.get_bromine().get_id()
                for fg in lig_structure.get_functional_groups()
            ]
            crest_data['binder_plane_deviations'] = [
                calculate_molecule_planarity(
                    mol=lig_structure.with_structure_from_file(i),
                    atom_ids=binder_ids,
                )
                for i in crest_conformer_files
            ]
            plot_binder_plane_deviation(
                measures=crest_data['binder_plane_deviations'],
                name=name,
                crest=True,
            )

        with open(crest_output_file, 'w') as f:
            json.dump(crest_data, f, indent=4)
        print('-------\n')


if __name__ == '__main__':
    main()
