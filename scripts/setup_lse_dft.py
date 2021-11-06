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

    charge = 0
    multiplicity = 1
    xyzfile = f'{base_name}_init.xyz'

    string = f'* xyzfile {charge} {multiplicity} {xyzfile}\n'

    struct.write(os.path.join(directory, xyzfile))

    return string


def write_input_file(input_file, dft_directory, mol_file):
    """
    Write ORCA single point energy of mol file.

    """

    top_line = (
        '! DFT SP RKS PBE0 def2-SVP D4 '
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

    xray_directory : (str)
        Directory with extracted-from-crystal ligands.
""")
        sys.exit()
    else:
        dft_directory = sys.argv[1]
        cage_directory = sys.argv[2]
        xray_directory = sys.argv[3]

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
    if os.path.exists(dft_directory):
        input(f'sure you want to delete {dft_directory}?')
        shutil.rmtree(dft_directory)
    os.mkdir(dft_directory)

    all_ligands = all_free_ligands+all_extr_ligands+all_xray_ligands
    all_input_file_names = []
    for lig_mol in sorted(all_ligands):
        if 'jd354' in lig_mol:
            continue
        input_file = lig_mol.replace('.mol', '.in').split('/')[-1]
        write_input_file(input_file, dft_directory, lig_mol)
        all_input_file_names.append(input_file)

    write_sub_file(all_input_file_names, dft_directory)
    print(len(all_input_file_names))


if __name__ == "__main__":
    main()
