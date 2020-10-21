#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to calculate flexibility measure of ligands.

Author: Andrew Tarzia

Date Created: 21 Oct 2020

"""

import sys
from os.path import exists, join
from glob import glob

import stk
import stko

from molecule_building import metal_FFs

from atools import (
    build_conformers,
    MOC_collapse_mc,
    MOC_uff_opt,
    calculate_molecule_planarity,
    update_from_rdkit_conf,
    get_atom_distance,
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
        print(f'doing {lig}...')
        lig_structure = ligands[lig]
        binder_plane_deviations = calculate_flex(lig_structure)
        plot_bpd(binder_plane_deviations, lig)

    sys.exit()


if __name__ == '__main__':
    main()
