#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to convert list of REFCODEs into CIFs - no constraints are applied.

Author: Andrew Tarzia

Date Created: 12 May 2019

"""
from ccdc.io import EntryReader, CrystalWriter
import sys


def main():
    if (not len(sys.argv) == 3):
        print """
    Usage: REFCODEs_to_CIFs.py REFCODE_file missing_struct
        REFCODE_file (str) - file with list of REFCODEs
        missing_struct (str) - file with list of REFCODEs with missing structs
        """
        sys.exit()
    else:
        RCODE_file = sys.argv[1]
        missing_struct = sys.argv[2]

    # read in CSD
    entry_reader = EntryReader('CSD')

    REFCODEs = []
    for line in open(RCODE_file, 'r'):
        REFCODEs.append(line.rstrip())

    count = 0
    count_no = 0
    RC_nostruct = []
    for i, RC in enumerate(sorted(REFCODEs)):
        count_no += 1
        entry = entry_reader.entry(RC)
        crystal = entry.crystal
        if entry.has_3d_structure is False:
            # print 'struct missing: '+str(RC)+' '+str(entry.ccdc_number)
            RC_nostruct.append(RC)
        else:
            # write to CIF
            CrystalWriter(RC+'_extracted.cif').write(crystal)
            count += 1
    print '-------------------------------------------------'
    print 'structures missing: '+str(len(RC_nostruct))+' of: '+str(len(REFCODEs))
    with open(missing_struct, 'w') as f:
        for RC in RC_nostruct:
            f.write(RC+'\n')


if __name__ == "__main__":
    main()
