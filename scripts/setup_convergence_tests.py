#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to setup parameter-search DFT calculations.

Author: Andrew Tarzia

Date Created: 20 Jan 2022

"""

import shutil
import sys
import glob
import os

import stk

import dft_utilities


def main():
    if (not len(sys.argv) == 3):
        print(
            """
Usage: setup_convergence_tests.py dft_directory cage_directory

    dft_directory : (str)
        Directory to run dft from - will be created.

    cage_directory : (str)
        Directory with required cage structures.

""")
        sys.exit()
    else:
        dft_directory = sys.argv[1]
        cage_directory = sys.argv[2]

    # Hardcode some test structures.
    all_cage_structures = [
        os.path.join(cage_directory, f'C_cl1_quad2_12_t_optc.mol'),
        os.path.join(cage_directory, f'C_cl1_quad2_3_c2v_optc.mol'),
        os.path.join(cage_directory, f'C_cl1_quad2_8_s62_optc.mol'),
    ]

    print(len(all_cage_structures))
    print(all_cage_structures)

    dft_directory = os.path.abspath(dft_directory)
    if os.path.exists(dft_directory):
        input(
            f'sure you want to delete {dft_directory}? ctrl-C if not!'
        )
        shutil.rmtree(dft_directory)
    os.mkdir(dft_directory)

    cutoffs = range(50, 501, 50)
    rel_cutoffs = range(10, 101, 10)
    for lig_mol in sorted(all_cage_structures):
        for cutoff in cutoffs:
            for rel_cutoff in rel_cutoffs:
                job_name = lig_mol.replace('.mol', '').split('/')[-1]
                job_name += f'_{cutoff}_{rel_cutoff}'
                spe = dft_utilities.CP2KEnergy(f'{job_name}_spe')
                molecule = stk.BuildingBlock.init_from_file(lig_mol)

                # Write single-point input file using coord from opt file.
                spe.write_calculation_input(
                    output_directory=dft_directory,
                    molecule=molecule,
                    charge=16,
                    guess='ATOMIC',
                    cutoff=cutoff,
                    rel_cutoff=rel_cutoff,
                    solvent=35.688,
                )
                # Write single-point slurm.
                spe.write_slurm_file(
                    output_directory=dft_directory,
                    hours=2,
                )


if __name__ == "__main__":
    main()
