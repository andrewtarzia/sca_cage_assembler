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
import glob
import os
sys.path.insert(0, '/home/atarzia/thesource/')
import pywindow_f
import IO_tools


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



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='')
    main()
