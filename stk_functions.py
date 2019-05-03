#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions that are useful for stk usage

Author: Andrew Tarzia

Date Created: 18 Mar 2019
"""

from glob import glob
import stk
import sys
from stk.molecular.molecules import MacroMoleculeBuildError
import json


def build_ABCBA(core, liga, link, flippedlink=False):
    '''Build ABCBA ligand using linear stk polymer.

    Keyword Arguments:
        core (stk.StructUnit) - molecule to use as core
        liga (stk.StructUnit) - molecule to use as liga
        link (stk.StructUnit) - molecule to use as link
        flippedlink (bool) - whether to flip the linker molecule (default False)

    Returns:
        polymer (stk.Polymer()) - polymer molecule pre optimization

    '''
    if flippedlink is False:
        orientation = [0, 0, 0, 1, 1]
    elif flippedlink is True:
        orientation = [0, 1, 0, 0, 1]
    polymer = stk.Polymer([liga, link, core],
                          stk.Linear(repeating_unit='ABCBA',
                                     orientation=orientation,
                                     n=1, ends='fg'))
    return polymer


def build_ABA(core, liga):
    '''Build ABCBA ligand using linear stk polymer.

    Keyword Arguments:
        core (stk.StructUnit) - molecule to use as core
        liga (stk.StructUnit) - molecule to use as liga
        link (stk.StructUnit) - molecule to use as link

    Returns:
        polymer (stk.Polymer()) - polymer molecule pre optimization

    '''
    polymer = stk.Polymer([liga, core],
                          stk.Linear(repeating_unit='ACA',
                                     orientation=[0, 0, 1],
                                     n=1, ends='fg'))
    return polymer


def build_population(directory, structunit, fgs=None, suffix='.mol'):
    '''Reads all SUFFIX files in directory into an stk population.

    Keyword Arguments:
        directory (str) - directory containing molecules
        structunit (str) - stk.StructUnit/2/3 class to use for molecule
        fgs (list) - list of functional groups on the molecules. Defaults to Br
        suffix (str) - file type (default .mol)

    Returns:
        popn (stk.Population()) - population of molecules

    '''
    if fgs is None:
        fgs = ['bromine']

    mols = []
    for file in sorted(glob(directory + '*' + suffix)):
        if structunit == 'StructUnit':
            mol = stk.StructUnit(file, fgs,
                                 name=file.rstrip(suffix).replace(directory, ''))
        elif structunit == 'StructUnit2':
            mol = stk.StructUnit2(file, fgs,
                                  name=file.rstrip(suffix).replace(directory, ''))
        elif structunit == 'StructUnit3':
            mol = stk.StructUnit3(file, fgs,
                                  name=file.rstrip(suffix).replace(directory, ''))
        mols.append(mol)
    popn = stk.Population(*mols)
    return popn


def topo_2_property(topology, property):
    '''Returns properties of a topology for a given topology name.

    Properties:
        'stk_func' - gives the stk topology function for building cages
        'stoich' - gives the stoichiometries of both building blocks assuming
            that the first building block has the larger number of functional groups.
        'noimines' - gives the number of imines formed to build that topology
        'expected_wind' - gives the number of windows expected

    Currently defined topologies:
        TwoPlusThree topologies
        ThreePlusThree topologies

    '''
    properties = ['stk_func', 'stoich', 'noimines', 'expected_wind']
    if property not in properties:
        print('{} not defined'.format(property))
        print('possible properties:', properties)
        sys.exit('exitting.')

    dict = {
        '2p3':
        {'stk_func': stk.two_plus_three.TwoPlusThree(),
         'stoich': (2, 3),
         'noimines': 6,
         'expected_wind': 3,
         },
        '4p6':
        {'stk_func': stk.two_plus_three.FourPlusSix(),
         'stoich': (4, 6),
         'noimines': 12,
         'expected_wind': 4,
         },
        '4p62':
        {'stk_func': stk.two_plus_three.FourPlusSix2(),
         'stoich': (4, 6),
         'noimines': 12,
         'expected_wind': 4,
         },
        '6p9':
        {'stk_func': stk.two_plus_three.SixPlusNine(),
         'stoich': (6, 9),
         'noimines': 18,
         'expected_wind': 5,
         },
        'dodec':
        {'stk_func': stk.two_plus_three.Dodecahedron(),
         'stoich': (20, 30),
         'noimines': 60,
         'expected_wind': 12,
         },
        '8p12':
        {'stk_func': stk.two_plus_three.EightPlusTwelve(),
         'stoich': (8, 12),
         'noimines': 24,
         'expected_wind': 6,
         },
        '1p1':
        {'stk_func': stk.three_plus_three.OnePlusOne(
         # place bb1 on vertex (0), bb2 on vertex (1)
         bb_positions={0: [0], 1: [1]}),
         'stoich': (1, 1),
         'noimines': 3,
         'expected_wind': 3,
         },
        '4p4':
        {'stk_func': stk.three_plus_three.FourPlusFour(
         # place bb1 on vertex (0, 2), bb2 on vertex (1, 3)
         bb_positions={0: [0, 3, 5, 6], 1: [1, 2, 4, 7]}),
         'stoich': (4, 4),
         'noimines': 12,
         'expected_wind': 6,
         },
    }
    if topology not in dict:
        print('properties not defined for {}'.format(topology))
        sys.exit('exitting')
    return dict[topology][property]


def is_collapse(topo, avg_diff, max_window_diam, cavity_size, no_window):
    expected_wind = topo_2_property(topo, property='expected_wind')
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


def load_StructUnitX(file, X=0):
    '''Load StructUnitX class with the cache turned off to avoid misreading of
    file.

    Keyword Arguments:
        file (str) - file to load in
        X (int) - order of struct unit to load. Defaults to 0

    X = 0 -> stk.StructUnit
    X = 2 -> stk.StructUnit2
    X = 3 -> stk.StructUnit3

    '''
    stk.OPTIONS['cache'] = False  # turn caching off for loading
    if X == 0:
        struct = stk.StructUnit(file)
    elif X == 2:
        struct = stk.StructUnit2(file)
    elif X == 3:
        struct = stk.StructUnit3(file)
    else:
        print('X must be 0, 2 or 3')
        sys.exit('exitting')
    stk.OPTIONS['cache'] = True  # turn caching back on
    return struct


def default_stk_MD_settings():
    '''Default settings from stk source code as of 26/04/19.

    '''
    Settings = {'output_dir': None,
                'timeout': None,
                'force_field': 16,
                'temperature': 300,  # K
                'conformers': 50,
                'time_step': 1.0,  # fs
                'eq_time': 10,  # ps
                'simulation_time': 200,  # ps
                'maximum_iterations': 2500,
                'minimum_gradient': 0.05,
                'use_cache': False
                }
    return Settings


def atarzia_short_MD_settings():
    '''My default settings for short, crude cage optimizations in stk.


    Modified on 26/04/19.
    '''
    Settings = {'output_dir': None,
                'timeout': None,
                'force_field': 16,
                'temperature': 700,  # K
                'conformers': 50,
                'time_step': 1,  # fs
                'eq_time': 50,  # ps
                'simulation_time': 1000,  # ps -- 1 ns
                'maximum_iterations': 2500,
                'minimum_gradient': 0.05,
                'use_cache': False}
    return Settings


def atarzia_long_MD_settings():
    '''My default settings for long and rigorous cage optimizations in stk.

    Mimics: Computationally-inspired discovery of an unsymmetrical porous organic cage
        DOI:10.1039/C8NR06868B

    Modified on 26/04/19.
    '''
    Settings = {'output_dir': None,
                'timeout': None,
                'force_field': 16,
                'temperature': 700,  # K
                'conformers': 10000,
                'time_step': 0.5,  # fs
                'eq_time': 500,  # ps
                'simulation_time': 100000,  # ps -- 100 ns
                'maximum_iterations': 2500,
                'minimum_gradient': 0.05,
                'use_cache': False}
    return Settings


def optimize_structunit(infile, outfile, exec,
                        settings=None, method='OPLS'):
    '''Read file into StructUnit and run optimization via method.

    '''
    # use standard settings applied in andrew_marsh work if md/settings is None
    if method == 'OPLS':
        if settings is None:
            Settings = default_stk_MD_settings()
        else:
            Settings = settings
        print(infile)
        struct = load_StructUnitX(infile, X=0)
        print('doing opt')
        # restricted=False optimization with OPLS forcefield by default
        ff = stk.MacroModelForceField(macromodel_path=exec, restricted=False)
        # MD process - run MD, collect N conformers, optimize each,
        # return lowest energy conformer
        md = stk.MacroModelMD(macromodel_path=exec,
                              output_dir=Settings['output_dir'],
                              timeout=Settings['timeout'],
                              force_field=Settings['force_field'],
                              temperature=Settings['temperature'],
                              conformers=Settings['conformers'],
                              time_step=Settings['time_step'],
                              eq_time=Settings['eq_time'],
                              simulation_time=Settings['simulation_time'],
                              maximum_iterations=Settings['maximum_iterations'],
                              minimum_gradient=Settings['minimum_gradient'],
                              use_cache=Settings['use_cache'])
        macromodel = stk.OptimizerSequence(ff, md)
        macromodel.optimize(mol=struct)
        struct.write(outfile)
        print(outfile)
    else:
        print('no other method is implemented yet.')
        sys.exit('exitting')


def get_OPLS3_energy_of_list(out_file, structures, macromod_,
                             dir='', opt=False, settings=None):
    '''Get OPLS3 single point energy of a list of structures.

    Keyword Arguments:
        file (str) - file that tracks the structures done to avoid redoing
        structures (list) - list of structure files
        dir (str) - directory where structures are if not in working dir
        opt (bool) - True if optimization is required using macromodel
        settings (dict) - settings for MacroModel Opt/MD
    '''
    # use standard settings applied in andrew_marsh work if md/settings is None
    if settings is None:
        Settings = default_stk_MD_settings()
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
            struct = load_StructUnitX(file, X=0)
            print('doing opt')
            # restricted=True optimization with OPLS forcefield by default
            ff = stk.MacroModelForceField(macromodel_path=exec,
                                          restricted=True)
            # MD process - run MD, collect N conformers, optimize each,
            # return lowest energy conformer
            md = stk.MacroModelMD(macromodel_path=exec,
                                  output_dir=Settings['output_dir'],
                                  timeout=Settings['timeout'],
                                  force_field=Settings['force_field'],
                                  temperature=Settings['temperature'],
                                  conformers=Settings['conformers'],
                                  time_step=Settings['time_step'],
                                  eq_time=Settings['eq_time'],
                                  simulation_time=Settings['simulation_time'],
                                  maximum_iterations=Settings['maximum_iterations'],
                                  minimum_gradient=Settings['minimum_gradient'],
                                  use_cache=Settings['use_cache'])
            macromodel = stk.OptimizerSequence(ff, md)
            macromodel.optimize(mol=struct)
            # get energy
            struct.energy.macromodel(16, macromod_)
            for i in struct.energy.values:
                energies[NAME] = struct.energy.values[i]
                calculated[NAME] = struct.energy.values[i]
        elif NAME not in calculated and opt is False:
            print('extracting energy')
            try:
                struct = load_StructUnitX(file, X=0)
            except TypeError:
                struct = load_StructUnitX(file + '_opt.mol', X=0)
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
                       settings=None, pdb=None):
    '''

    Keyword Arguments:
        prefix (str) - output file name prefix
        BB1 (str) - name of building block 1 file
        BB2 (str) - name of building block 2 file
        topology (stk.topology) = cage toplogy object
        macromod_ (str) - location of macromodel
        settings (dict) - settings for MacroModel Opt
        xyz (bool) - otuput XYZ file of optimized cage (default is None)
    '''
    # use standard settings applied in andrew_marsh work if md/settings is None
    if settings is None:
        Settings = default_stk_MD_settings()
    else:
        Settings = settings
    # try:
    cage = stk.Cage([BB1, BB2], topology)
    cage.write(prefix + '.mol')
    cage.dump(prefix + '.json')
    # restricted=True optimization with OPLS forcefield by default
    ff = stk.MacroModelForceField(macromodel_path=macromod_,
                                  restricted=True)
    # MD process - run MD, collect N conformers, optimize each,
    # return lowest energy conformer
    # no restricted
    md = stk.MacroModelMD(macromodel_path=macromod_,
                          output_dir=Settings['output_dir'],
                          timeout=Settings['timeout'],
                          force_field=Settings['force_field'],
                          temperature=Settings['temperature'],
                          conformers=Settings['conformers'],
                          time_step=Settings['time_step'],
                          eq_time=Settings['eq_time'],
                          simulation_time=Settings['simulation_time'],
                          maximum_iterations=Settings['maximum_iterations'],
                          minimum_gradient=Settings['minimum_gradient'],
                          use_cache=Settings['use_cache'])
    macromodel = stk.OptimizerSequence(ff, md)
    macromodel.optimize(mol=cage)
    cage.write(prefix + '_opt.mol')
    cage.dump(prefix + '_opt.json')
    if pdb is True:
        cage.write(prefix + '_opt.pdb')
    return cage
    # except MacroMoleculeBuildError:
    #     print('build failed')
    #     pass
