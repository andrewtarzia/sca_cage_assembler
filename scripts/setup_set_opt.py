#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to setup DFT optimisations of a cage set.

Author: Andrew Tarzia

Date Created: 15 Feb 2021

"""

import shutil
import sys
import glob
import os

import stk

import dft_utilities


def main():
    if (not len(sys.argv) == 5):
        print(
            """
Usage: setup_lse_dft.py dft_directory cage_directory xray_directory
delete_directory

    dft_directory : (str)
        Directory to run dft from - will be created.

    cage_directory : (str)
        Directory with required cage structures.

    cage_set : (str)
        Name of cage set to use in DFT calculations.

    delete_directory : (str)
        't' to delete dft_directory, other to not.

""")
        sys.exit()
    else:
        dft_directory = sys.argv[1]
        cage_directory = sys.argv[2]
        cage_set = sys.argv[3]
        delete_directory = True if sys.argv[4] == 't' else False

    all_cage_structures = glob.glob(os.path.join(
        cage_directory, f'C_{cage_set}_*optc.mol'
    ))

    print(len(all_cage_structures))

    dft_directory = os.path.abspath(dft_directory)
    if os.path.exists(dft_directory) and delete_directory:
        input(
            f'sure you want to delete {dft_directory}? ctrl-C if not!'
        )
        shutil.rmtree(dft_directory)
    if not os.path.exists(dft_directory):
        os.mkdir(dft_directory)

    chosen_cutoff = 700
    chosen_rel_cutoff = 60
    for lig_mol in sorted(all_cage_structures):
        job_name = lig_mol.replace('.mol', '').split('/')[-1]
        opt = dft_utilities.CP2KOptimizer(f'{job_name}_opt')
        molecule = stk.BuildingBlock.init_from_file(lig_mol)

        # Write optimisation input file.
        opt.write_calculation_input(
            output_directory=dft_directory,
            molecule=molecule,
            charge=16,
            guess='ATOMIC',
            cutoff=chosen_cutoff,
            rel_cutoff=chosen_rel_cutoff,
        )
        # Write optimisation slurm.
        opt.write_slurm_file(
            output_directory=dft_directory,
            hours=12,
        )


if __name__ == "__main__":
    main()
