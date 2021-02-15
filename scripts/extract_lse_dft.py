#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to extract results from DFT calculations of ligand strain.

Author: Andrew Tarzia

Date Created: 15 Feb 2021

"""

import sys
import os
import glob
import pandas as pd

from atools import read_gfnx2xtb_eyfile


def read_orca_output(output_path):

    with open(output_path, 'r') as f:
        lines = f.readlines()

    targ_str = 'FINAL SINGLE POINT ENERGY'
    for line in lines:
        if targ_str in line:
            return float(line.rstrip().split(' ')[-1])


def collect_all_energies(ey_files, dft_directory):

    TEST_FILES = [
        'C_cl1_quad2_16_th2_sg114_1_1.ey',
        'jd235_sg120_1_0.ey',
        'jd235_sg120_1_5.ey'
    ]

    energies = {}
    for xtb_ey_file in ey_files:
        test = xtb_ey_file.split('/')[-1]
        if test not in TEST_FILES:
            continue
        xtb_ey = read_gfnx2xtb_eyfile(xtb_ey_file)
        output_file = xtb_ey_file.replace('.ey', '.out').split('/')[-1]
        print(xtb_ey_file, xtb_ey, output_file)
        output_path = os.path.join(dft_directory, output_file)
        if not os.path.exists(output_path):
            raise FileNotFoundError(xtb_ey_file)
        dft_ey = read_orca_output(output_path)
        if dft_ey is None:
            raise ValueError(output_path)
        print(xtb_ey, dft_ey)
        energies[xtb_ey_file] = (output_file, xtb_ey, dft_ey)

    return energies


def calculate_all_strain_energies(all_energies):

    strain_energies = {}

    for file in all_energies:
        if 'opt' in file:
            continue
        extr_xtb_energy = all_energies[file][1]
        extr_dft_energy = all_energies[file][2]

        opt_file = 'something'
        free_xtb_energy = all_energies[opt_file][1]
        free_dft_energy = all_energies[opt_file][2]

        xtb_strain = extr_xtb_energy - free_xtb_energy
        dft_strain = extr_dft_energy - free_dft_energy
        if 'jd' in file:
            from_ = 'crystal'
        else:
            from_ = 'computation'
        strain_energies[file] = (xtb_strain, dft_strain, from_)

    return strain_energies


def main():
    if (not len(sys.argv) == 4):
        print(
            """Usage: extract_lse_dft.py dft_directory cage_directory
                            xray_directory

    dft_directory : (str)
        Directory to run dft from - will be created.

    cage_directory : (str)
        Directory with required cage structures.

    xray_directory : (str)
        Directory with extracted ligands."""
        )
        sys.exit()
    else:
        dft_directory = sys.argv[1]
        cage_directory = sys.argv[2]
        xray_directory = sys.argv[3]

    xtb_cage_ey_files = [
        i for i in glob.glob(os.path.join(cage_directory, '*.ey'))
        if 'optc.ey' not in i
    ]
    xtb_xray_ey_files = glob.glob(os.path.join(xray_directory, '*.ey'))

    print(len(xtb_cage_ey_files))
    print(len(xtb_xray_ey_files))

    all_energies = collect_all_energies(
        ey_files=xtb_cage_ey_files+xtb_xray_ey_files,
        dft_directory=dft_directory,
    )
    print(all_energies)

    strain_energies = calculate_all_strain_energies(all_energies)
    print(strain_energies)

    df = pd.DataFrame.from_dict({
        'file': [i for i in strain_energies],
        'xtb': [strain_energies[i][0] for i in strain_energies],
        'dft': [strain_energies[i][1] for i in strain_energies],
        'from': [strain_energies[i][2] for i in strain_energies],
    })
    df.to_csv('strain_energy_comparison.csv')


if __name__ == "__main__":
    main()
