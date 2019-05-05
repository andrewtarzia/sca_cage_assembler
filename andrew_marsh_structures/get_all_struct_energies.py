#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run MD and GFN optimizations of cages or precursors to get energies.

Author: Andrew Tarzia

Date Created: 03 May 2019
"""

import sys
from os.path import join, isfile
from glob import glob
import pandas as pd
sys.path.insert(0, '/home/atarzia/thesource/')
from GFN_functions import run_GFN_base, get_energies
from stk_functions import optimize_structunit, atarzia_long_MD_settings
from IO_tools import convert_MOL3000_2_PDB_XYZ


def main():
    """Run script.

    """
    if (not len(sys.argv) == 4):
        print("""
Usage: get_all_struct_energies.py output_file suffix
    output_file: file to output results
    suffix (str) - file suffix to run GFN calculation on - should end in .mol
    MD (str) - 't' if a macromodel conformer search should be applied
""")
        sys.exit()
    else:
        output_file = sys.argv[1]
        suffix = sys.argv[2]
        MD = sys.argv[3]

    macromod_ = '/home/atarzia/software/schrodinger_install'

    if isfile(output_file):
        # read in data
        DATA = pd.read_csv(output_file)
    else:
        DATA = pd.DataFrame(columns=['files', 'NAMES', 'FE (au)'])

    done_files = list(DATA['files'])
    # iterate over all files with suffix
    if suffix != '_opt.mol':
        # not on cages
        files = [i for i in glob('*' + suffix) if '_opt.mol' not in i]
    else:
        # being run on cages
        files = [i for i in glob('*' + suffix)]
    # ignore done files
    files = [i for i in files
             if i.replace('.mol', '_opt.mol') not in done_files]
    print('files to do GFN:', files)
    print('files done:', done_files)
    if MD.lower() == 't':
        # do a macromodelMD conformer search with stk first
        outfiles = []
        for file in files:
            OUTFILE = file.replace('.mol', '_opt.mol')
            if isfile(OUTFILE) is False:
                optimize_structunit(infile=file,
                                    outfile=OUTFILE,
                                    exec=macromod_,
                                    settings=atarzia_long_MD_settings())
            outfiles.append(OUTFILE)
    else:
        outfiles = [i.replace('.mol', '_opt.mol') for i in files
                    if i.replace('.mol', '_opt.mol') not in done_files]

    print('output files to analyze:', outfiles)
    # do GFN optimization
    # no need to rerun structures already in DATA
    xyzs = [i.replace('.mol', '.xyz') for i in outfiles
            if i.rstrip('.mol') not in list(DATA['NAMES'])]
    print('XYZs to do:', xyzs)
    # convert outfiles (.mol) to XYZ files (.xyz)
    for file in outfiles:
        if file.rstrip('.xyz') not in list(DATA['NAMES']):
            convert_MOL3000_2_PDB_XYZ(file=file)

    failed = run_GFN_base(xyzs=xyzs)
    if len(failed) > 0:
        sys.exit('---> Some GFN calcs failed. Exitting.')

    # get Free energies and save to output_file
    for file in outfiles:
        DIR = file.replace('.mol', '')
        OUT = file.replace('.mol', '.output')
        energy_results = get_energies(output_file=join(DIR, OUT),
                                      GFN_exec='/home/atarzia/software/xtb_190418/bin/xtb')
        DATA = DATA.append({'NAMES': file.rstrip('.mol'),
                            'files': file,
                            'FE (au)': energy_results['FE']},
                           ignore_index=True)
        # write dataFrame
        DATA.to_csv(output_file, index=False)


if __name__ == "__main__":
    main()
