#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Minimum working example of bug with latest stk version.

Author: Andrew Tarzia

Date Created: 08 Apr 2019

"""

import sys
import stk


def main():
    if (not len(sys.argv) == 1):
        print("""
    Usage: bug_mwe.py
        """)
        sys.exit()
    else:
        # CIF = sys.argv[1]
        pass

    # load in molecules
    core = stk.StructUnit('core.mol', ['bromine'])
    liga = stk.StructUnit('liga.mol', ['bromine'])
    link = stk.StructUnit('link.mol', ['bromine'])

    # build ABCBA molecule
    polymer = stk.Polymer([liga, link, core],
                          stk.Linear(repeating_unit='ABCBA',
                                     orientation=[0, 0, 0, 1, 1],
                                     n=1, ends='fg'))
    # output as built
    json_file = 'Test.json'
    polymer.dump(json_file)
    mol_file = 'Test.mol'
    polymer.write(mol_file)


if __name__ == "__main__":
    main()
