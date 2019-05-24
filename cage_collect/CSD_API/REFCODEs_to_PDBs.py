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
import sys
sys.path.insert(0, '/home/atarzia/thesource/')
import CSD_f


def main():
    if (not len(sys.argv) == 4):
        print """
    Usage: REFCODEs_to_CIFs.py REFCODE_file missing_struct
        REFCODE_file (str) - file with list of REFCODEs
        missing_struct (str) - file with list of REFCODEs with missing structs
        cross_references (str) - file with list of REFCODEs that require cross_references
        """
        sys.exit()
    else:
        RCODE_file = sys.argv[1]
        missing_struct = sys.argv[2]
        cross_references = sys.argv[3]

    # read in CSD and updates
    entry_reader = CSD_f.get_entryreader()

    REFCODEs = []
    for line in open(RCODE_file, 'r'):
        REFCODEs.append(line.rstrip())

    RC_nostruct = []
    RC_CR = []
    for i, RC in enumerate(sorted(REFCODEs)):
        print 'doing: '+str(RC)
        entry = entry_reader.entry(RC)
        crystal = None
        if entry.has_3d_structure:
            crystal = entry.crystal
        elif entry.has_3d_structure is False:
            # test if CSD REFCODE is of type XXXXXX01
            # which implies that XXXXXX will have coordinates and this is a
            # child entry
            # only assuming this can be the case a new REFCODE is in
            if len(entry.cross_references) == 0:
                # print 'struct missing: '+str(RC)+' '+str(entry.ccdc_number)
                RC_nostruct.append(RC)
                continue
            else:
                for CR in entry.cross_references:
                    # check if cross ref type is coordinates
                    if CR.type == 'Coordinates ref':
                        idents = CR.identifiers
                        for ID in idents:
                            try:
                                new_entry = entry_reader.entry(ID)
                            except RuntimeError:
                                # implies this new entry ID is not in the CSD
                                RC_nostruct.append(RC)
                                continue
                            if new_entry.has_3d_structure:
                                crystal = new_entry.crystal
                                RC_CR.append((RC, ID))
                                break
        # write to CIF - saves as REFCODE in input file even if cross reference
        # is used
        if crystal is not None:
            packed = crystal.packing()
            ccdc.io.CrystalWriter(RC+'_extracted.pdb').write(packed)
            # use ASE to add cell parameters and resave PDB
            CELL = [crystal.cell_lengths.a, crystal.cell_lengths.b,
                    crystal.cell_lengths.c, crystal.cell_angles.alpha,
                    crystal.cell_angles.beta, crystal.cell_angles.gamma]
            CSD_f.rewrite_pdb(pdb=RC+'_extracted.pdb', cell=CELL)
    print '-------------------------------------------------'
    print 'structures missing: '+str(len(RC_nostruct))+' of '+str(len(REFCODEs))
    with open(missing_struct, 'w') as f:
        for RC in RC_nostruct:
            f.write(RC+'\n')
    print '-------------------------------------------------'
    print 'cross refs used: '+str(len(RC_CR))+' of '+str(len(REFCODEs))
    with open(cross_references, 'w') as f:
        for RC, CR in RC_CR:
            f.write(RC+','+CR+'\n')


if __name__ == "__main__":
    main()
