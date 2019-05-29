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
        """
        sys.exit()
    else:
        file_suffix = sys.argv[1]



if __name__ == "__main__":
    main()
