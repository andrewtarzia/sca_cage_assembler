#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to calculate flexibility measure of ligands.

Author: Andrew Tarzia

Date Created: 21 Oct 2020

"""

import numpy as np
import sys
from os.path import join, exists
from glob import glob
import json

import stk

from atools import (
    build_conformers,
    calculate_molecule_planarity,
    colors_i_like,
    histogram_plot_N,
    update_from_rdkit_conf,
)


def load_ligands(directory):

    ligands = {}
    for lig in glob(join(directory, f'quad2*_opt.mol')):
        l_name = lig.replace(directory, '').replace('_opt.mol', '')
        ligands[l_name] = stk.BuildingBlock.init_from_file(
            lig,
            functional_groups=[stk.BromoFactory()],
        )

    return ligands


def calculate_flex(molecule):
    """
    Calculate flexibility as plane deviation of binder groups.

    Assumes molecule has N functional groups with only 1 binder atom
    each.

    """

    for fg in molecule.get_functional_groups():
        if len(list(fg.get_bonder_ids())) > 1:
            raise ValueError(
                f'{molecule} has functional groups with more'
                ' than 1 binder.'
            )

    bromo_ids = [
        fg.get_bromine().get_id()
        for fg in molecule.get_functional_groups()
    ]

    binder_plane_deviations = []
    cids, confs = build_conformers(
        mol=molecule,
        N=200,
        ETKDG_version='v3'
    )
    new_molecule = molecule.clone()
    for cid in cids:
        # Update stk_mol to conformer geometry.
        new_molecule = update_from_rdkit_conf(
            stk_mol=new_molecule,
            rdk_mol=confs,
            conf_id=cid
        )

        binder_plane_deviations.append(
            calculate_molecule_planarity(
                mol=new_molecule,
                atom_ids=bromo_ids,
            )
        )

    return binder_plane_deviations


def plot_bpd(measures, name):

    fig, ax = histogram_plot_N(
        Y=measures,
        X_range=(0, 30),
        width=0.2,
        alpha=1.0,
        color=colors_i_like()[1],
        edgecolor=colors_i_like()[1],
        xtitle=r'binder atom plane deviation [$\mathrm{\AA}$]',
        N=1
    )
    fig.tight_layout()
    fig.savefig(
        f'{name}_planedev_dist.pdf',
        dpi=720,
        bbox_inches='tight'
    )


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

    for lig in sorted(ligands):
        lig_structure = ligands[lig]
        json_file = f'{lig}_planedev_dist.json'
        print(f'doing {lig}...')
        if exists(json_file):
            with open(json_file, 'r') as f:
                binder_plane_deviations = json.load(f)
        else:
            binder_plane_deviations = calculate_flex(lig_structure)
            with open(json_file, 'w') as f:
                json.dump(binder_plane_deviations, f)
        print(
            lig,
            max(binder_plane_deviations),
            np.std(binder_plane_deviations)
        )
        plot_bpd(binder_plane_deviations, lig)

    sys.exit()


if __name__ == '__main__':
    main()
