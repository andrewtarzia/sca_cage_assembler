#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to sort CIFs or PDBs based on their pyWindow results (<- extract_indep_cages.py)

Author: Andrew Tarzia

Date Created: 05 Apr 2019

"""

import logging
import sys
import pandas as pd
import glob
import os
sys.path.insert(0, '/home/atarzia/thesource/')
import pywindow_f
import IO_tools


def main():
    if (not len(sys.argv) == 4):
        print("""
    Usage: sort_structures.py DB_file output_file file_type
        DB_file (str) - file with initial list of REFCODEs
        output_file (str) - file to output results of sorting to
        file_type (str) - set whether to run on PDBs (enter: pdb) or CIFs (enter: cif)
        """)
        sys.exit()
    else:
        DB_file = sys.argv[1]
        output_file = sys.argv[2]
        file_type = sys.argv[3]

    # temporary check for non-implemented issue with extractedm.cif cases
    # these cases were manually collected
    for i in glob.glob(f'*.{file_type}'):
        if 'extractedm' in i:
            logging.error('This code cannot handled extractedm cases! Please implement them.')

    refcodes = sorted([i.rstrip() for i in open(DB_file, 'r').readlines()])
    files = [i+'_extracted.'+file_type for i in refcodes]
    logging.info(f'> started with: {len(refcodes)} structures to sort.')
    if os.path.isfile(output_file):
        # read structures already checked to avoid double calculations
        OUTDATA = pd.read_csv(output_file)
        done_files = list(OUTDATA.file)
        logging.info(f'> {len(done_files)} structures already done.')
    else:
        # write output file
        with open(output_file, 'w') as f:
            f.write('file,deleted\n')
        OUTDATA = pd.read_csv(output_file)
        done_files = []

    # iterate over files
    count = len(done_files)
    for file in files:
        # skip done structures
        if file in done_files:
            continue
        if os.path.isfile(file):
            if file_type == 'cif':
                pdb = IO_tools.convert_CIF_2_PDB(file, wstruct=False)
            elif file_type == 'pdb':
                pdb = IO_tools.check_ASE_handle(file, wstruct=False)
            if pdb is None:
                logging.warning(f'> ASE failed to load {file}')
                # file failed to load in ASE
                OUTDATA = OUTDATA.append({'file': file, 'deleted': 'M'},
                                         ignore_index=True)
                os.remove(file)
            else:
                logging.info(f'> doing {count} of {len(files)}')
                # check if at least one molecule has a pore_diameter_opt > 0.25 angstrom
                if pywindow_f.check_PDB_for_pore(file=pdb, diam=0.0):
                    OUTDATA = OUTDATA.append({'file': file, 'deleted': 'N'},
                                             ignore_index=True)
                else:
                    # delete molecule if not
                    OUTDATA = OUTDATA.append({'file': file, 'deleted': 'Y'},
                                             ignore_index=True)
                    os.remove(file)
                    try:
                        os.remove(pdb)
                    except FileNotFoundError:
                        pass
                    os.remove(pdb.replace('.pdb', '_rebuild.pdb'))
        else:
            # file missing.
            OUTDATA = OUTDATA.append({'file': file, 'deleted': 'M'},
                                     ignore_index=True)

        # add to done cifs
        done_files.append(file)
        # update output file
        OUTDATA.to_csv(output_file, index=False)
        count += 1

    remaining = list(OUTDATA[OUTDATA['deleted'] == 'N']['file'])
    logging.info(f'> ended with: {len(remaining)} structures.')


if __name__ == "__main__":
    # logging.basicConfig(level=logging.DEBUG, format='%(levelname)s-%(message)s')
    # logging.debug(f'Debug mode!')
    logging.basicConfig(level=logging.INFO, format='')
    main()
