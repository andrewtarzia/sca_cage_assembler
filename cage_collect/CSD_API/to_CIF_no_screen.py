#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to convert list of REFCODEs into CIFs without screening.

Author: Andrew Tarzia

Date Created: 4 Mar 2019

"""
from ccdc.io import EntryReader, CrystalWriter

# read in CSD
entry_reader = EntryReader('CSD')

RCODE_file = 'CSD_search_FEBUpdate.gcd'
print('reading', RCODE_file, 'is that correct??')
REFCODEs = []
for line in open(RCODE_file, 'r'):
    REFCODEs.append(line.rstrip())

count = 0
count_no = 0
RC_list = []
for i, RC in enumerate(sorted(REFCODEs)):
    print RC
    # if RC.lower() != 'dovmuw':
    #     continue
    count_no += 1
    entry = entry_reader.entry(RC)
    # print('%s %s' % (RC, entry.ccdc_number))
    # print(entry.chemical_name)
    # print(entry.is_polymeric)
    # skip polymeric structures
    if entry.chemical_name is not None:
        if 'catena' in entry.chemical_name:
            continue
    if entry.is_polymeric is True:
        continue
    # skip if structure is powder study
    if entry.is_powder_study is True:
        continue
    # break
    crystal = entry.crystal
    if entry.has_3d_structure is False:
        print RC, entry.ccdc_number
        RC_list.append(RC)
    else:
        # write to CIF
        CrystalWriter(RC+'_extracted.cif').write(crystal)
        count += 1

print(count, 'cifs found from', count_no, 'RCs')
print RC_list
