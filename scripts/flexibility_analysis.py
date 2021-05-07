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
from os.path import join, exists
from glob import glob
import matplotlib.pyplot as plt
import json

import stk

from atools import (
    calculate_molecule_planarity,
    colors_i_like,
    histogram_plot_N,
    crest_conformer_search,
    get_atom_distance,
)

from utilities import split_xyz_file


def load_ligands(directory):

    ligands = {}
    for lig in glob(join(directory, 'quad2*_opt.mol')):
        l_name = lig.replace(directory, '').replace('_opt.mol', '')
        ligands[l_name] = stk.BuildingBlock.init_from_file(
            lig,
            functional_groups=[stk.BromoFactory()],
        )

    return ligands


def calculate_flex(molecule, name, la_pairs):
    """
    Calculate flexibility of molecule from CREST conformer ensemble.

    Three approaches:
        1) Based on plane deviation of Bromine atoms in molecule.
        Assumes molecule has N functional groups with only 1 binder
        atom each.

        2) Based on plane deviation of all atoms in molecule.

        3) Long-axis distance.

    """

    for fg in molecule.get_functional_groups():
        if len(list(fg.get_bonder_ids())) > 1:
            raise ValueError(
                f'{molecule} has functional groups with more'
                ' than 1 binder.'
            )

    # Crest part.
    if not exists(f'crst_{name}/crest_rotamers.xyz'):
        new_molecule = crest_conformer_search(
            molecule=molecule,
            output_dir=f'crst_{name}',
            gfn_exec='/home/atarzia/anaconda3/envs/sca_cages/bin/xtb',
            crest_exec='/home/atarzia/software/crest/crest',
            gfn_version=2,
            nc=3,
            opt_level='crude',
            charge=0,
            keepdir=False,
            cross=False,
            etemp=300,
            no_unpaired_e=0,
            speed_setting='squick',
            solvent=('acetonitrile', 'normal'),
        )
        new_molecule.write(f'{name}_loweconf.mol')

    # Extract some measure of conformer ensemble size.
    crest_output_file = f'{name}_flex_measure.json'
    crest_data = {}
    with open(f'crst_{name}/crest.output', 'r') as f:
        for line in f.readlines():
            # Get number of conformers.
            if ' number of unique conformers for further calc' in line:
                crest_data['no_conformers'] = (
                    int(line.rstrip().split(' ')[-1])
                )

            # Get number of rotamers.
            if 'total number unique points considered further' in line:
                crest_data['no_rotamers'] = (
                    int(line.rstrip().split(' ')[-1])
                )

    # Analyse all rotamers from CREST.
    crest_conformer_files = split_xyz_file(
        num_atoms=molecule.get_num_atoms(),
        xyz_file=f'crst_{name}/crest_conformers.xyz',
    )
    bromo_ids = [
        fg.get_bromine().get_id()
        for fg in molecule.get_functional_groups()
    ]
    crest_data['binder_plane_deviations'] = [
        calculate_molecule_planarity(
            mol=molecule.with_structure_from_file(i),
            atom_ids=bromo_ids,
        )
        for i in crest_conformer_files
    ]
    crest_data['plane_deviations'] = [
        calculate_molecule_planarity(
            mol=molecule.with_structure_from_file(i),
        )
        for i in crest_conformer_files
    ]
    print(f'{name} has {len(crest_conformer_files)} conformers')

    pair1_cents = [
        molecule.with_structure_from_file(i).get_centroid(
            atom_ids=tuple(la_pairs[0])
        )
        for i in crest_conformer_files
    ]
    pair2_cents = [
        molecule.with_structure_from_file(i).get_centroid(
            atom_ids=tuple(la_pairs[1])
        )
        for i in crest_conformer_files
    ]

    la_dist = [
        np.linalg.norm(i-j)
        for i, j in zip(pair1_cents, pair2_cents)
    ]
    crest_data['la_dist'] = la_dist

    with open(crest_output_file, 'w') as f:
        json.dump(crest_data, f)

    return crest_data


def plot_lap(measures, name, crest=False):
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


def plot_pd(measures, name, crest=False):

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


def plot_bpd(measures, name, crest=False):

    fig, ax = histogram_plot_N(
        Y=measures,
        X_range=(0, 30),
        width=0.2,
        alpha=1.0,
        color=colors_i_like()[1],
        edgecolor=colors_i_like()[1],
        xtitle=r'all atom plane deviation [$\mathrm{\AA}$]',
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


def get_long_axis_atoms(ligands):

    long_axis_atom_pairs = {}
    for ligand in ligands:
        molecule = ligands[ligand]
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

        pairs = (
            [i[0] for i in short_vector_fg_1],
            [i[0] for i in short_vector_fg_2]
        )
        long_axis_atom_pairs[ligand] = pairs

    return long_axis_atom_pairs


def main():
    first_line = (
        'Usage: flexibility_analysis.py '
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

    # Load in each ligand structure.
    ligands = load_ligands(ligand_directory)
    long_axis_atom_pairs = get_long_axis_atoms(ligands)

    for lig in sorted(ligands):
        lig_structure = ligands[lig]
        crest_data = calculate_flex(
            molecule=lig_structure,
            name=lig,
            la_pairs=long_axis_atom_pairs[lig],
        )

        print(
            '::',
            abs(max(crest_data['la_dist'])-min(crest_data['la_dist']))
        )
        plot_lap(crest_data['la_dist'], lig, crest=True)
        plot_bpd(
            crest_data['binder_plane_deviations'],
            lig,
            crest=True
        )
        plot_pd(crest_data['plane_deviations'], lig, crest=True)
        print(
            lig,
            crest_data['no_rotamers'],
            crest_data['no_conformers']
        )
        print('---')


if __name__ == '__main__':
    main()
