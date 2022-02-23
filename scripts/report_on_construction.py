#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to report on successful and failed constructions/optimisations.

Author: Andrew Tarzia

Date Created: 23 Feb 2022

"""

import sys
import json
import os
import glob
import stk


def main():
    first_line = (
        'Usage: report_on_construction.py'
    )
    if (not len(sys.argv) == 1):
        print(f"""
{first_line}
    """)
        sys.exit()
    else:
        pass

    cs_files = glob.glob('*CS.json')
    _output_directory = 'construction_report'
    _output_string = ''

    if not os.path.exists(_output_directory):
        os.mkdir(_output_directory)

    count_optc_files = 0
    count_optimized = 0
    count_unoptimized = 0
    count_expected = 0
    for cage_set_file in cs_files:
        cage_set = cage_set_file.replace('_CS.json', '')
        _output_string += '-----------------\n'
        _output_string += (
            f'{cage_set} from {cage_set_file} report:\n'
        )
        with open(cage_set_file, 'r') as f:
            cs_data = json.load(f)

        for cage_name in cs_data:
            count_expected += 1
            cage_data = cs_data[cage_name]
            symm = cage_name.split('_')[-1]
            _output_string += f'\nSymmetry: {symm}\n'
            opt_file = f'{cage_name}_optc.mol'
            if os.path.exists(opt_file):
                _output_string += f'1) {opt_file} exists\n'
                count_optc_files += 1
                stk.BuildingBlock.init_from_file(opt_file).write(
                    os.path.join(
                        _output_directory,
                        opt_file.replace('.mol', '.pdb')
                    )
                )
            else:
                _output_string += f'1) {opt_file} does not exist\n'

            if cage_data['optimized']:
                count_optimized += 1
                _output_string += f'2) optimized!\n'
            else:
                count_unoptimized += 1
                _output_string += f'2) not optimized!\n'

        _output_string  += '-----------------\n'

    _output_string += (
        'In summary:\n'
        f'Count expected: {count_expected}\n'
        f'Count optc files: {count_optc_files}\n'
        f'Count optimized: {count_optimized}\n'
        f'Count unoptimized: {count_unoptimized}\n'
        f'Expected met?: '
        f'{count_optimized+count_unoptimized == count_expected}\n'
    )

    with open(os.path.join(_output_directory, 'report.txt'), 'w') as f:
        f.write(_output_string)

    print(_output_string)


if __name__ == "__main__":
    main()
