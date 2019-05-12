#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to sort CIFs based on their pyWindow results (<- extract_indep_cages.py)

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
    if (not len(sys.argv) == 3):
        print("""
    Usage: sort_structures.py DB_file output_file
        DB_file (str) - file with initial list of REFCODEs
        output_file (str) - file to output results of sorting to
        """)
        sys.exit()
    else:
        DB_file = sys.argv[1]
        output_file = sys.argv[2]

    # temporary check for non-implemented issue with extractedm.cif cases
    # these cases were manually collected
    for i in glob.glob('*.cif'):
        if 'extractedm' in i:
            logging.error('This code cannot handled extractedm cases! Please implement them.')

    refcodes = sorted([i.rstrip() for i in open(DB_file, 'r').readlines()])
    cifs = [i+'_extracted.cif' for i in refcodes]
    logging.info(f'> started with: {len(refcodes)} CIFs to sort.')
    if os.path.isfile(output_file):
        # read CIFs already checked to avoid double calculations
        OUTDATA = pd.read_csv(output_file)
        done_cifs = list(OUTDATA.cif)
        logging.info(f'> {len(done_cifs)} CIFs already done.')
    else:
        # write output file
        with open(output_file, 'w') as f:
            f.write('cif,deleted\n')
        OUTDATA = pd.read_csv(output_file)
        done_cifs = []

    # iterate over CIFs
    count = len(done_cifs)
    for cif in cifs:
        # skip done cifs
        if cif in done_cifs:
            continue
        if os.path.isfile(cif):
            pdb = IO_tools.convert_CIF_2_PDB(cif, wstruct=False)
            logging.info(f'> doing {count} of {len(cifs)}')
            # check if at least one molecule has a pore_diameter_opt > 0.25 angstrom
            if pywindow_f.check_PDB_for_pore(file=pdb, diam=0.25):
                OUTDATA = OUTDATA.append({'cif': cif, 'deleted': 'N'},
                                         ignore_index=True)
            else:
                # delete molecule if not
                OUTDATA = OUTDATA.append({'cif': cif, 'deleted': 'Y'},
                                         ignore_index=True)
                os.remove(cif)
                os.remove(pdb)
                os.remove(pdb.replace('.pdb', '_rebuild.pdb'))
        else:
            # CIF missing.
            OUTDATA = OUTDATA.append({'cif': cif, 'deleted': 'M'},
                                     ignore_index=True)

        # add to done cifs
        done_cifs.append(cif)
        # update output file
        OUTDATA.to_csv(output_file, index=False)
        count += 1

    remaining = list(OUTDATA[OUTDATA['deleted'] == 'N'])
    logging.info(f'> ended with: {len(remaining)} CIFs.')


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, format='%(levelname)s-%(message)s')
    logging.debug(f'Debug mode!')
    # logging.basicConfig(level=logging.INFO, format='')
    main()
