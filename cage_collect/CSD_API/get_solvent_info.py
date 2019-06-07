#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to convert list of REFCODEs into PDBs. No symm and constraints are applied.

Author: Andrew Tarzia

Date Created: 24 May 2019

"""
import ccdc.io
from ase.io import read
import glob
import sys
sys.path.insert(0, '/home/atarzia/thesource/')
import CSD_f


def main():
    if (not len(sys.argv) == 2):
        print """
    Usage: get_solvent_info.py file_suffix
        file_suffix (str) - suffix following REFCODE of files to analyse
            ('_extracted.pdb')
        """
        sys.exit()
    else:
        file_suffix = sys.argv[1]

    # read in CSD and updates
    entry_reader = CSD_f.get_entryreader()

    files = sorted(glob.glob('*'+file_suffix))
    REFCODEs = []
    for file in files:
        REFCODEs.append(file.replace(file_suffix, ''))

    RC_properties = []
    for i, RC in enumerate(sorted(REFCODEs)):
        print 'doing: '+str(RC)
        entry = entry_reader.entry(RC)
        # note structures with solvent
        # solvent = 'n'
        # if entry.chemical_name is not None:
        #     if len(entry.chemical_name.split(' ')) > 1:
        #         solvent = 'y'
        # disorder details
        squeeze = 'n'
        disorder = 'n'
        DD = entry.disorder_details
        if DD is not None:
            disorder = 'y'
            if 'squeeze' in DD.lower():
                squeeze = 'y'
            elif 'platon' in DD.lower():
                squeeze = 'y'
            elif 'mask' in DD.lower():
                squeeze = 'y'
        print DD
        rem = entry.remarks
        print rem
        RC_prop = (RC, squeeze, disorder)
        print entry.chemical_name
        print 'squeeze:'+squeeze
        print 'disorder:'+disorder
        print RC_prop
        # sys.exit()
        RC_properties.append(RC_prop)
        # a = raw_input('')

    # write to file
    with open('solvent_prop.csv', 'w') as f:
        for i in RC_properties:
            f.write(i[0]+','+i[1]+','+i[2]+'\n')
    print '-------------------------------------------------'


if __name__ == "__main__":
    main()
