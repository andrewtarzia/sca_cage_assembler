#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions that are useful for CSD API.

Author: Andrew Tarzia

Date Created: 12 May 2019
"""

import ccdc.io
from ase.io import read
import glob
import os


def get_entryreader():
    '''Get entry reader and updates.

    Example from
    https://downloads.ccdc.cam.ac.uk/documentation/API/modules/io_api.html#module-ccdc.io
    - .sqlite is known file format, .inf is not.

    '''
    directory = ccdc.io.csd_directory()
    csd_and_updates = glob.glob(os.path.join(directory, '*.sqlite'))
    csd_and_updates_reader = ccdc.io.EntryReader(csd_and_updates)
    return csd_and_updates_reader


def rewrite_pdb(pdb, cell):
    '''Overwrite structure in pdb with cell parameters.

    '''
    s = read(pdb)
    s.set_pbc(True)
    a, b, c, alpha, beta, gamma = cell
    s.set_cell([a, b, c, alpha, beta, gamma])
    # wrap atoms
    s.wrap()
    s.write(pdb)
