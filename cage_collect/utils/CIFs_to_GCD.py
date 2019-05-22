#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to take all REFCODEs of CIFs in a dir and write a .gcd for ConQuest.

Author: Andrew Tarzia

Date Created: 22 May 2019

"""

import logging
import sys
import glob


def main():
    if (not len(sys.argv) == 2):
        print("""
    Usage: CIFs_to_GCD.py gcd_file
        gcd_file (str) - file to output REFCODEs to
        """)
        sys.exit()
    else:
        gcd_file = sys.argv[1]

    # temporary check for non-implemented issue with extractedm.cif cases
    # these cases were manually collected
    for i in glob.glob('*.cif'):
        if 'extractedm' in i:
            logging.error('This code cannot handled extractedm cases! Please implement them.')

    CIFs = sorted(glob.glob('*.cif'))
    output_str = ''
    for cif in CIFs:
        RC = cif.replace('_extracted.cif', '')
        output_str += RC+'\n'

    with open(gcd_file, 'w') as f:
        f.write(output_str)
    logging.info(f'> wrote gcd with: {len(CIFs)} CIFs.')


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='')
    main()
