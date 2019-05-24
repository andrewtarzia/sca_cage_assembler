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
import os
sys.path.insert(0, '/home/atarzia/thesource/')
import pywindow_f


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
    pdbs = [i+'_extracted.pdb' for i in refcodes]
    logging.info(f'> started with: {len(refcodes)} structures to classify.')

    if os.path.isfile(output_file):
        # read CIFs already checked to avoid double calculations
        OUTDATA = pd.read_csv(output_file)
        done_RCs = list(OUTDATA['REFCODE'])
        logging.info(f'> {len(done_RCs)} structures already done.')
    else:
        # write output file
        with open(output_file, 'w') as f:
            f.write('REFCODE,molecule,pore_diam_opt,no_windows\n')
        OUTDATA = pd.read_csv(output_file)
        done_RCs = []

    # iterate over CIFs
    count = len(done_RCs)
    for pdb in pdbs:
        RC = pdb.replace('_extracted.pdb', '')
        # skip done cifs
        if RC in done_RCs:
            continue
        if os.path.isfile(pdb) is False:
            raise(f'{pdb} not present!')
        logging.info(f'> doing {count} of {len(pdbs)}: {RC}')
        # load and modularize pdb
        rbs = pywindow_f.modularize(file=pdb)
        if rbs is None:
            # handle pyWindow failure
            raise(f'{pdb} failed modularize!')
        # iterate over all molecules, skipping those with n_atoms < 5
        Mol = rbs.molecules
        for molec in Mol:
            mol = Mol[molec]
            if mol.no_of_atoms < 5:
                continue
            # run analysis
            try:
                analysis = mol.full_analysis()
            except ValueError:
                logging.warning(f'{pdb}_{molec} failed pywindow full_analysis.')
                analysis = None
            # define output
            if analysis is None:
                continue
            pdo = analysis['pore_diameter_opt']['diameter']
            if analysis['windows']['diameters'] is not None:
                nwind = len(analysis['windows']['diameters'])
            else:
                nwind = 0
            # if it is a cage:
            if pdo > 0.0 and nwind >= 2:
                # add to output
                OUTDATA = OUTDATA.append({'REFCODE': RC, 'molecule': molec,
                                          'pore_diam_opt': pdo,
                                          'no_windows': nwind},
                                         ignore_index=True)
                # output structure
                Mol[molec].dump_molecule(
                    RC + "_MP_{0}_coms.pdb".format(molec),
                    include_coms=True,
                    override=True)
                Mol[molec].dump_molecule(
                    RC + "_MP_{0}.pdb".format(molec),
                    include_coms=False,
                    override=True)
        # add to done cifs
        done_RCs.append(RC)
        # update output file
        OUTDATA.to_csv(output_file, index=False)
        count += 1


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='')
    main()
