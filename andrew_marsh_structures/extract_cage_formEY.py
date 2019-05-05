#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to calculate the formation energies of a list of cages.

Author: Andrew Tarzia

Date Created: 03 May 2019
"""

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
if __name__ == "__main__":
    main()
