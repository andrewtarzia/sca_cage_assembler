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

    if isfile(output_file):
        # read in data
        DATA = pd.read_csv(output_file)
    else:
        DATA = pd.DataFrame(columns=['files', 'NAMES', 'FE (au)'])
    # get Free energies and save to output_file
    for file in xyzs:
        DIR = file.replace('.xyz', '')
        OUT = file.replace('.xyz', '.out')
        ORIG_FILE = file.replace('.xyz', '.mol')
        energy_results = get_energies(output_file=join(DIR, OUT))
        DATA = DATA.append({'NAMES': file.rstrip('.xyz'),
                            'files': ORIG_FILE,
                            'FE (au)': energy_results['FE']},
                           ignore_index=True)
        # write dataFrame
        DATA.to_csv(output_file, index=False)


if __name__ == "__main__":
    main()
