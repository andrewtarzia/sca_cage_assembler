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


def write_sub_file(input_files, dft_directory):
    orca_bin_dir = '/apps/orca/4.2.1/bin/orca'

    runfile = os.path.join(f'{dft_directory}', 'run_orca.sh')

    runlines = ''.join([
        f"{orca_bin_dir} {infile} > {infile.replace('.in', '.out')}\n"
        for infile in input_files
    ])

    string = (
        '#PBS -N lig_spe\n'
        '#PBS -l walltime=72:00:00\n'
        '#PBS -l select=1:ncpus=32:mem=124gb\n\n'
        'module load orca/4.2.1\n\n'
        'cd $PBS_O_WORKDIR\n\n'
        f'{runlines}'
    )

    with open(runfile, 'w') as f:
        f.write(string)


def write_molecule_section(directory, base_name, struct):

    charge = 16
    multiplicity = 1
    xyzfile = f'{base_name}_init.xyz'

    string = f'* xyzfile {charge} {multiplicity} {xyzfile}\n'
    string += (
        f'#* xyzfile {charge} {multiplicity} {base_name}_B97-3c\n'
    )

    struct.write(os.path.join(directory, xyzfile))

    return string


def write_input_file(input_file, dft_directory, mol_file):
    """
    Write ORCA optimisation of mol file.

    """

    top_line = (
        '! DFT COPT B97-3c '
        'Grid6 NOFINALGRID SlowConv TightSCF\n\n'
    )
    base_name = '_'+input_file.replace('.in', '')
    base_line = f'%base "{base_name}" \n\n'

    scf_section = (
        '%scf\n   MaxIter 2000\nend\n\n'
    )

    procs_section = (
        '%pal\n   nprocs 32\nend\n\n'
    )

    struct = stk.BuildingBlock.init_from_file(mol_file)
    mol_section = write_molecule_section(
        dft_directory, base_name, struct
    )

    string = top_line
    string += base_line
    string += scf_section
    string += procs_section
    string += mol_section

    with open(f'{dft_directory}/{input_file}', 'w') as f:
        f.write(string)


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
        input_file = lig_mol.replace('.mol', '.in').split('/')[-1]
        write_input_file(input_file, dft_directory, lig_mol)
        all_input_file_names.append(input_file)

    write_sub_file(all_input_file_names, dft_directory)
    print(len(all_input_file_names))


if __name__ == "__main__":
    main()
