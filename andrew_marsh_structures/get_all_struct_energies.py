#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to compare the properties of cages from the same amine but two different
aldehydes.

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
    files = glob('*' + suffix)
    # ignore done files
    files = [i for i in files if i not in done_files]
    if MD.lower() == 't':
        # do a macromodelMD conformer search with stk first
        outfiles = []
        for file in files:
            OUTFILE = file.replace('.mol', '_opt.mol')
            optimize_structunit(infile=file,
                                outfile=OUTFILE,
                                exec=macromod_,
                                settings=atarzia_long_MD_settings())
            outfiles.append(OUTFILE)
    else:
        outfiles = files

    # do GFN optimization
    xyzs = [i.replace('.mol', '.xyz') for i in outfiles]
    failed = run_GFN_base(xyzs=xyzs)
    if len(failed) > 0:
        sys.exit('---> Some GFN calcs failed. Exitting.')

    # get Free energies and save to output_file
    for file in xyzs:
        DIR = file.replace('.xyz', '')
        OUT = file.replace('.xyz', '.out')
        ORIG_FILE = file.replace('.xyz', '.mol')
        energy_results = get_energies(output_file=join(DIR, OUT))
        DATA = DATA.append({'NAMES': file.rstrip('.xyz'),
                            'files': ORIG_FILE,
                            'FE (au)': energy_results['FE']},
                           ignore_index=True)
        # write dataFrame
        DATA.to_csv(output_file, index=False)


if __name__ == "__main__":
    main()
