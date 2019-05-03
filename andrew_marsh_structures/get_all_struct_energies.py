#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to compare the properties of cages from the same amine but two different
aldehydes.

Author: Andrew Tarzia

Date Created: 01 May 2019
"""

import sys
sys.path.insert(0, '/home/atarzia/thesource/')
def main():
    """Run script.

    """
    if (not len(sys.argv) == 4):
        print("""
Usage: get_all_struct_energies.py output_file suffix
    output_file: file to output results
    suffix (str) - file suffix to run GFN calculation on - should end in .mol
    MD (str) - 't' if a macromodel conformer search should be applied
""")
        sys.exit()
    else:
        output_file = sys.argv[1]
        suffix = sys.argv[2]
        MD = sys.argv[3]

    macromod_ = '/home/atarzia/software/schrodinger_install'

if __name__ == "__main__":
    main()
