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

    # iterate over CIFs
    count = len(done_cifs)
    for cif in cifs:
        # skip done cifs
        if cif in done_cifs:
            continue
        if os.path.isfile(cif):
            pdb = IO_tools.convert_CIF_2_PDB(cif, wstruct=False)
            if pdb is None:
                logging.warning(f'> ASE failed to load {cif}')
                OUTDATA = OUTDATA.append({'cif': cif, 'BUG?': 'M'},
                                         ignore_index=True)
            else:
                logging.info(f'> doing {count} of {len(cifs)}')
                # check if at least one molecule has a pore_diameter_opt > 0.25 angstrom
                if has_bug(pdb):
                    OUTDATA = OUTDATA.append({'cif': cif, 'BUG?': 'Y'},
                                             ignore_index=True)
                else:
                    # delete molecule if not
                    OUTDATA = OUTDATA.append({'cif': cif, 'BUG?': 'N'},
                                             ignore_index=True)

        # add to done cifs
        done_cifs.append(cif)
        # update output file
        OUTDATA.to_csv(output_file, index=False)
        count += 1

    wbug = list(OUTDATA[OUTDATA['BUG?'] == 'Y']['cif'])
    wASEbug = list(OUTDATA[OUTDATA['BUG?'] == 'M']['cif'])
    logging.info(f'> ended with: {len(wbug)} buggy CIFs.')
    logging.info(f'> ended with: {len(wASEbug)} buggy CIFs with ASE.')


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='')
    main()
