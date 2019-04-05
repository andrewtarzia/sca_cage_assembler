#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to sort CIFs based on their pyWindow results (<- extract_indep_cages.py)

Author: Andrew Tarzia

Date Created: 05 Apr 2019

"""

import sys
from glob import glob
from os.path import isfile
import json
sys.path.insert(0, '/home/atarzia/thesource/')


def main():
    all_cifs = sorted(glob('*extracted.cif')) + sorted(glob('*extractedm.cif'))

    for cif in all_cifs:
        cif_out_1 = cif.replace('.cif', '_cleaned1.cif')
        cif_out_2 = cif.replace('.cif', '_cleaned_P1.cif')
        RC = cif.split('_')[0]
        print(RC)
        RC_files = glob(RC + '_*')
        print(RC_files)
        # remove those that got skipped in structure_preparation.py
        if isfile(cif_out_1) is False and isfile(cif_out_2) is False:
            continue
        # remove those without .json files
        json_files = glob(RC + '*.json')
        if not len(json_files):
            continue
        for js in json_files:
            print(js)
            with open(js, 'r') as f:
                data = json.load(f)
            print(data['pore_diameter_opt']['diameter'])
            print(data['pore_volume_opt'])
            print(data['windows']['diameters'])
        break


if __name__ == "__main__":
    main()
