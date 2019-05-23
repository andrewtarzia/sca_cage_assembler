#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to extract the most porous molecule from a rebuilt pdb.

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
    Usage: extract_most_porous.py DB_file output_file
        DB_file (str) - file with initial list of REFCODEs
        output_file (str) - file to output results of sorting to
        """)
        sys.exit()
    else:
        DB_file = sys.argv[1]
        output_file = sys.argv[2]



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='')
    main()
