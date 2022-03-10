#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to setup DFT calculations of cage stability.

Author: Andrew Tarzia

Date Created: 10 Mar 2022

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

    dft_directory = os.path.abspath(dft_directory)

    chosen_cutoff = 700
    chosen_rel_cutoff = 60
    for lig_mol in sorted(all_cage_structures):
        job_name = lig_mol.replace('.mol', '').split('/')[-1]
        spe = dft_utilities.CP2KEnergy(f'{job_name}_spe')
        raise SystemExit('Load structure from optimised coordinates.')
        raise SystemExit('Provide error if output from opt not there.')
        molecule = stk.BuildingBlock.init_from_file(lig_mol)

        # Write single-point input file using coord from opt file.
        spe.write_calculation_input(
            output_directory=dft_directory,
            molecule=molecule,
            charge=16,
            guess='ATOMIC',
            cutoff=chosen_cutoff,
            rel_cutoff=chosen_rel_cutoff,
        )
        # Write single-point slurm.
        spe.write_slurm_file(
            output_directory=dft_directory,
            hours=2,
        )
        raise SystemExit()


if __name__ == "__main__":
    main()
