#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to assist in the classification of CIFs based on the molecules in the CIF.

Author: Andrew Tarzia

Date Created: 23 May 2019

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
    Usage: classify_structures.py DB_file output_file
        DB_file (str) - file with initial list of REFCODEs
        output_file (str) - file to output results of sorting to
        """)
        sys.exit()
    else:
        DB_file = sys.argv[1]
        output_file = sys.argv[2]

    refcodes = sorted([i.rstrip() for i in open(DB_file, 'r').readlines()])
    rebuilt_pdbs = [i+'_extracted_rebuild.pdb' for i in refcodes]
    logging.info(f'> started with: {len(refcodes)} structures to classify.')
    if os.path.isfile(output_file):
        # read CIFs already checked to avoid double calculations
        OUTDATA = pd.read_csv(output_file)
        done_RCs = list(OUTDATA['REFCODE'])
        logging.info(f'> {len(done_RCs)} structures already done.')
    else:
        # write output file
        with open(output_file, 'w') as f:
            f.write('REFCODE,molecule,pore_diam_opt,no_windows,classification\n')
        OUTDATA = pd.read_csv(output_file)
        done_RCs = []

    # iterate over CIFs
    count = len(done_RCs)
    for rbpdb in rebuilt_pdbs:
        RC = rbpdb.replace('_extracted_rebuild.pdb', '')
        # skip done cifs
        if RC in done_RCs:
            continue
        if os.path.isfile(rbpdb):
            logging.info(f'> doing {count} of {len(RC)}: {RC}')
            # load in rebuilt structure

            # modularize

            # iterate over all molecules, skipping those with n_atoms < 5
            for molec in MOLECULES:
                # run analysis

                # define output
                pdo = 0
                nwind = 0

                # output structure

                # visualize in PROGRAM to get class

                CL = 'M'

                # close program


                OUTDATA = OUTDATA.append({'REFCODE': RC, 'molecule': molec,
                                          'pore_diam_opt': pdo,
                                          'no_windows': nwind,
                                          'classification': CL},
                                         ignore_index=True)

        else:
            # pdb missing.
            molec = 0
            pdo = 0
            nwind = 0
            CL = 'M'  # mistake
            OUTDATA = OUTDATA.append({'REFCODE': RC, 'molecule': molec,
                                      'pore_diam_opt': pdo,
                                      'no_windows': nwind,
                                      'classification': CL},
                                     ignore_index=True)

        # add to done cifs
        done_RCs.append(RC)
        # update output file
        OUTDATA.to_csv(output_file, index=False)
        count += 1


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='')
    main()
