#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run GFN jobs for all XYZ files in a directory.

Author: Andrew Tarzia

Date Created: 08 Mar 2019

"""

import sys
from glob import glob
sys.path.insert(0, '/home/atarzia/thesource/')
import GFN_f


def main():
    if (not len(sys.argv) == 2):
        print("""
Usage: run_GFN.py suffix
    suffix (str) - file suffix to run GFN calculation on - should end in .xyz
    """)
        sys.exit()
    else:
        suffix = sys.argv[1]
    xyzs = glob('*' + suffix)
    GFN_f.GFN_from_xyzs(xyzs=xyzs)


if __name__ == "__main__":
    main()
