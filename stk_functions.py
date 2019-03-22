#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions that are useful for stk usage

Author: Andrew Tarzia

Date Created: 18 Mar 2019
"""

import stk
from stk.molecular.molecules import MacroMoleculeBuildError
import json


def expected_window(topo):
    e_wind = {'dodec': 12, '4p6': 4, '4p62': 4,
              '8p12': 6,  '6p9': 5, '2p3': 3,
              '4p4': 6, '1p1': 3, '2p2': 4}
    return e_wind[topo]


def is_collapse(topo, avg_diff, max_window_diam, cavity_size, no_window):
    expected_wind = expected_window(topo)
    if expected_wind == no_window:
        alpha = 4 * avg_diff / (max_window_diam * expected_wind)
        if alpha < 0.035 and cavity_size > 1:
            # not collapsed
            return False
        else:
            # unknown
            return None
    else:
        # collapsed
        return True


def get_asymmetry(data):
    """Calculate assymetry as defined in GA paper (Berardo)

    The sum of all the windows' pair differences represents the asymmetry
    of the individual, Asymmetry parameter in eqn (1)

    Deprecated 22/03/19 -- Andrew Tarzia
    """
    print('you should not be using this.')
    window_sizes = data['windows']['diameters']
    total = 0
    for i, a in enumerate(window_sizes):
        for j, b in enumerate(window_sizes[i:]):
            if i != j+i:
                diff = abs(a - b)
                total += diff
    return total


def get_OPLS3_energy_of_list(out_file, structures, macromod_,
                             dir='', opt=False, md=None, settings=None):
    '''Get OPLS3 single point energy of a list of structures.

    Keyword Arguments:
        file (str) - file that tracks the structures done to avoid redoing
        structures (list) - list of structure files
        dir (str) - directory where structures are if not in working dir
        opt (bool) - True if optimization is required using macromodel
        md (dict) - settings for MacroModel MD
        settings (dict) - settings for MacroModel Opt
    '''
    # use standard settings applied in andrew_marsh work if md/settings is None
    if md is None:
        MD = {'timeout': None,
              'force_field': 16,
              'temp': 700,
              'confs': 50,
              'time_step': 1.0,
              'eq_time': 10,
              'sim_time': 200,
              'max_iter': 2500,
              'gradient': 0.05}
    else:
        MD = md
    if settings is None:
        Settings = {'restricted': False,
                    'timeout': None,
                    'force_field': 16,
                    'max_iter': 2500,
                    'gradient': 0.05,
                    'md': True}
    else:
        Settings = settings
    with open(out_file, 'r') as f:
        calculated = json.load(f)
    energies = {}
    for file in structures:
        print(file)
        NAME = file.replace(dir, '').replace('.mol', '')
        # optimize
        if NAME not in calculated and opt is True:
            struct = stk.StructUnit(file)
            print('doing opt')
            stk.macromodel_opt(struct,
                               macromodel_path=macromod_,
                               settings=Settings,
                               md=MD)
            # get energy
            struct.energy.macromodel(16, macromod_)
            for i in struct.energy.values:
                energies[NAME] = struct.energy.values[i]
                calculated[NAME] = struct.energy.values[i]
        elif NAME not in calculated and opt is False:
            print('extracting energy')
            try:
                struct = stk.StructUnit(file)
            except TypeError:
                struct = stk.StructUnit(file+'_opt.mol')
            struct.energy.macromodel(16, macromod_)
            for i in struct.energy.values:
                energies[NAME] = struct.energy.values[i]
                calculated[NAME] = struct.energy.values[i]
        else:
            print('already calculated')
            energies[NAME] = calculated[NAME]
    # save energies of calculated precursors to avoid recalculation
    with open(out_file, 'w') as f:
        json.dump(calculated, f)
    return energies


def build_and_opt_cage(prefix, BB1, BB2, topology, macromod_,
                       md=None, settings=None):
    '''

    Keyword Arguments:
        prefix (str) - output file name prefix
        BB1 (str) - name of building block 1 file
        BB2 (str) - name of building block 2 file
        topology (stk.topology) = cage toplogy object
        macromod_ (str) - location of macromodel
        md (dict) - settings for MacroModel MD
        settings (dict) - settings for MacroModel Opt
    '''
    # use standard settings applied in andrew_marsh work if md/settings is None
    if md is None:
        MD = {'timeout': None,
              'force_field': 16,
              'temp': 700,
              'confs': 50,
              'time_step': 1.0,
              'eq_time': 10,
              'sim_time': 200,
              'max_iter': 2500,
              'gradient': 0.05}
    else:
        MD = md
    if settings is None:
        Settings = {'restricted': False,
                    'timeout': None,
                    'force_field': 16,
                    'max_iter': 2500,
                    'gradient': 0.05,
                    'md': True}
    try:
        cage = stk.Cage([BB1, BB2], topology)
        cage.write(prefix+'.mol')
        stk.macromodel_opt(cage, macromodel_path=macromod_,
                           settings=Settings,
                           md=MD)
        cage.write(prefix+'_opt.mol')
        return cage
    except MacroMoleculeBuildError:
        pass
