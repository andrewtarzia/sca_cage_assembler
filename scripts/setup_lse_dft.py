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

    xray_directory : (str)
        Directory with extracted-from-crystal ligands.

    delete_directory : (str)
        't' to delete dft_directory, other to not.

""")
        sys.exit()
    else:
        dft_directory = sys.argv[1]
        cage_directory = sys.argv[2]
        xray_directory = sys.argv[3]
        delete_directory = True if sys.argv[4] == 't' else False

    all_free_ligands = glob.glob(os.path.join(
        cage_directory, '*_sg*_opt.mol'
    ))
    all_extr_ligands = [
        i for i in glob.glob(os.path.join(
            cage_directory, '*_sg*.mol'
        ))
        if i not in all_free_ligands
    ]
    all_xray_ligands = glob.glob(os.path.join(
        xray_directory, '*_sg*.mol'
    ))

    print(len(all_free_ligands))
    print(len(all_extr_ligands))
    print(len(all_xray_ligands))

    dft_directory = os.path.abspath(dft_directory)
    if os.path.exists(dft_directory) and delete_directory:
        input(
            f'sure you want to delete {dft_directory}? ctrl-C if not!'
        )
        shutil.rmtree(dft_directory)
    if not os.path.exists(dft_directory):
        os.mkdir(dft_directory)

    all_ligands = all_free_ligands+all_extr_ligands+all_xray_ligands
    all_input_file_names = []
    for lig_mol in sorted(all_ligands):
        input_file = lig_mol.replace('.mol', '.in').split('/')[-1]
        dft_utilities.write_input_file(
            input_file=input_file,
            dft_directory=dft_directory,
            mol_file=lig_mol,
            charge=0,
            method='spe-pbe'
        )
        all_input_file_names.append(input_file)

    dft_utilities.write_sub_file(all_input_file_names, dft_directory)
    print(len(all_input_file_names))


if __name__ == "__main__":
    main()
