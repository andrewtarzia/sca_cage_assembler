#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to setup DFT calculations of ligand strain.

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
    if (not len(sys.argv) == 4):
        print(
            """
Usage: setup_lse_dft.py dft_directory cage_directory xray_directory

    dft_directory : (str)
        Directory to run dft from - will be created.

    cage_directory : (str)
        Directory with required cage structures.

    cage_set : (str)
        Name of cage set to use in DFT calculations.

""")
        sys.exit()
    else:
        dft_directory = sys.argv[1]
        cage_directory = sys.argv[2]
        cage_set = sys.argv[3]

    all_cage_structures = glob.glob(os.path.join(
        cage_directory, f'C_{cage_set}_*optc.mol'
    ))

    print(len(all_cage_structures))
    print(all_cage_structures)

    dft_directory = os.path.abspath(dft_directory)
    if os.path.exists(dft_directory):
        input(
            f'sure you want to delete {dft_directory}? ctrl-C if not!'
        )
        shutil.rmtree(dft_directory)
    os.mkdir(dft_directory)

    all_input_file_names = []
    for lig_mol in sorted(all_cage_structures):
        job_name = lig_mol.replace('.mol', '').split('/')[-1]
        opt = dft_utilities.CP2KOptimizer(f'{job_name}_opt')
        spe = dft_utilities.CP2KEnergy(f'{job_name}_spe')
        molecule = stk.BuildingBlock.init_from_file(lig_mol)

        # Write optimisation input file.
        opt.write_calculation_input(
            output_directory=dft_directory,
            molecule=molecule,
            charge=16,
            guess='ATOMIC',
            cutoff=350,
            rel_cutoff=60,
            solvent=35.688,
        )
        # Write single-point input file using coord from opt file.
        spe.write_calculation_input(
            output_directory=dft_directory,
            molecule=None,
            charge=16,
            guess='ATOMIC',
            cutoff=350,
            rel_cutoff=60,
            solvent=35.688,
        )
        # Write optimisation slurm.
        opt.write_slurm_file(
            output_directory=dft_directory,
            hours=12,
        )
        # Write single-point slurm.
        spe.write_slurm_file(
            output_directory=dft_directory,
            hours=2,
        )
        raise SystemExit()

        dft_utilities.write_input_file(
            input_file=input_file,
            dft_directory=dft_directory,
            mol_file=lig_mol,
            charge=16,
            method='opt-b97',
        )
        all_input_file_names.append(input_file)

    dft_utilities.write_sub_file(all_input_file_names, dft_directory)
    print(len(all_input_file_names))


if __name__ == "__main__":
    main()
