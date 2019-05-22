#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to test for the AZOBEY bug in all CIFs from POC_list_2

Author: Andrew Tarzia

Date Created: 22 May 2019

"""
def main():
    if (not len(sys.argv) == 3):
        print("""
    Usage: test_azobeybug.py DB_file output_file
        DB_file (str) - file with initial list of REFCODEs
        output_file (str) - file to output results of sorting to
        """)
        sys.exit()
    else:
        DB_file = sys.argv[1]
        output_file = sys.argv[2]

    refcodes = sorted([i.rstrip() for i in open(DB_file, 'r').readlines()])
    cifs = [i+'_extracted.cif' for i in refcodes]
    if os.path.isfile(output_file):
        # read CIFs already checked to avoid double calculations
        OUTDATA = pd.read_csv(output_file)
        done_cifs = list(OUTDATA.cif)
        logging.info(f'> {len(done_cifs)} CIFs already done.')
    else:
        # write output file
        with open(output_file, 'w') as f:
            f.write('cif,BUG?\n')
        OUTDATA = pd.read_csv(output_file)
        done_cifs = []

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='')
    main()
