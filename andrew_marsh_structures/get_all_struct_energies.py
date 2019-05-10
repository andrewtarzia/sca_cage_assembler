#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run MD and GFN optimizations of cages or precursors to get energies.

Author: Andrew Tarzia

Date Created: 03 May 2019
"""

import sys
import os
import glob
import pandas as pd
sys.path.insert(0, '/home/atarzia/thesource/')
import GFN_f
import stk_f
import IO_tools


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

    if os.path.isfile(output_file):
        # read in data
        DATA = pd.read_csv(output_file)
    else:
        DATA = pd.DataFrame(columns=['files', 'NAMES', 'FE (au)'])

    done_files = list(DATA['files'])
    # iterate over all files with suffix
    if suffix != '_opt.mol':
        # not on cages
        files = [i for i in glob.glob('*' + suffix) if '_opt.mol' not in i]
        # ignore done files
        files = [i for i in files
                 if i.replace('.mol', '_opt.mol') not in done_files]
    else:
        # being run on cages
        files = [i for i in glob.glob('*' + suffix)]
        # ignore done files
        files = [i for i in files
                 if i not in done_files]
    print('files to do GFN:', files)
    print('files done:', done_files)
    if MD.lower() == 't':
        # do a macromodelMD conformer search with stk first
        outfiles = []
        for file in files:
            OUTFILE = file.replace('.mol', '_opt.mol')
            if os.path.isfile(OUTFILE) is False:
                stk_f.optimize_structunit(infile=file,
                                                  outfile=OUTFILE,
                                                  exec=macromod_,
                                                  settings=stk_f.atarzia_long_MD_settings())
            outfiles.append(OUTFILE)
    else:
        if suffix != '_opt.mol':
            outfiles = [i.replace('.mol', '_opt.mol') for i in files
                        if i.replace('.mol', '_opt.mol') not in done_files]
        else:
            outfiles = [i for i in files
                        if i not in done_files]

    print('output files to analyze:', outfiles)
    # do GFN optimization
    # no need to rerun structures already in DATA
    xyzs = [i.replace('.mol', '.xyz') for i in outfiles
            if i.rstrip('.mol') not in list(DATA['NAMES'])]
    print('XYZs to do:', xyzs)
    # convert outfiles (.mol) to XYZ files (.xyz)
    for file in outfiles:
        if file.rstrip('.xyz') not in list(DATA['NAMES']):
            IO_tools.convert_MOL3000_2_PDB_XYZ(file=file)

    failed = GFN_f.run_GFN_base(xyzs=xyzs)
    if len(failed) > 0:
        print('--------------------------------------------------------')
        print('---> Some GFN calcs failed.')
        print('--------------------------------------------------------')

    # get Free energies and save to output_file
    for file in xyzs:
        if file in failed:
            continue
        DIR = file.replace('.xyz', '')
        OUT = file.replace('.xyz', '.output')
        energy_results = GFN_f.get_energies(output_file=os.path.join(DIR, OUT),
                                            GFN_exec='/home/atarzia/software/xtb_190418/bin/xtb')
        if 'FE' in energy_results:
            DATA = DATA.append({'NAMES': file.rstrip('.xyz'),
                                'files': file,
                                'FE (au)': float(energy_results['FE'])},
                                ignore_index=True)
        # write dataFrame
        DATA.to_csv(output_file, index=False)


if __name__ == "__main__":
    main()
