#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run a particular script in batch on a series of files in current DIR.

Author: Andrew Tarzia

Date Created: 24 Apr 2019

"""

import sys
from glob import glob


def main():
    if (not len(sys.argv) == 4):
        print("""
Usage: batch_processor.py script glob_pattern warn
    script (str) - full path to python script to run
    glob_pattern (str) - pattern to use for glob
    warn (str) - if warn is 't/T' the program continues over failed outputs
        will raise error otherwise
    """)
        sys.exit()
    else:
        script = sys.argv[1]
        glob_pattern = sys.argv[2]
        warn = sys.argv[3].lower()

    command_items = ['python', script]

    files = glob(glob_pattern)
    print('{} files to run on'.format(len(files)))
    count = 1
    for file in files:
        print('doing {} of {}'.format(count, len(files)))
        command_items += file
        print(command_items)
        # execute and determine the result



        print('done')


if __name__ == "__main__":
    main()
