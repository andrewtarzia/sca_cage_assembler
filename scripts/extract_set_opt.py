#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to extract structures from DFT optimisations of a cage set.

Author: Andrew Tarzia

Date Created: 22 Feb 2022

"""

import shutil
import sys
import glob
import os

import stk

import dft_utilities


def split_xyz_file(num_atoms, xyz_file):
    """
    Splits xyz trajectory file into xyz files.

    """

    with open(xyz_file, 'r') as f:
        lines = f.readlines()

    file_strings = []
    string = []
    for line in lines:
        if f' {num_atoms} ' in f' {line.strip()} ':
            if len(string) == 0:
                string.append(line)
            else:
                # New block.
                file_strings.append(string)
                string = [line]
        else:
            string.append(line)
    # Add last set.
    file_strings.append(string)

    return file_strings


def main():
    if (not len(sys.argv) == 4):
        print(
            """
Usage: extract_set_opt.py dft_directory cage_directory xray_directory
delete_directory

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

    for lig_mol in sorted(all_cage_structures):
        opt_file = lig_mol.split('/')[-1]
        opt_mol = stk.BuildingBlock.init_from_file(opt_file)
        opt_xyz_traj = os.path.join(
            dft_directory,
            opt_file.replace('.mol', '_opt-pos-1.xyz')
        )
        opt_xyz_file = os.path.join(
            dft_directory,
            opt_file.replace('.mol', '_opt-pos-1_fin.xyz')
        )

        traj_xyz = split_xyz_file(
            num_atoms=opt_mol.get_num_atoms(),
            xyz_file=opt_xyz_traj,
        )

        final_step = traj_xyz[-1]
        final_energy = float(final_step[1].split()[-1])
        print(f'{opt_file} final energy: {final_energy} a.u.')

        with open(opt_xyz_file, 'w') as f:
            for line in final_step:
                f.write(line)

        new_opt_file = opt_file.replace('optc', 'optdft')
        new_opt_mol = opt_mol.with_structure_from_file(opt_xyz_file)
        new_opt_mol.write(new_opt_file)


if __name__ == "__main__":
    main()
