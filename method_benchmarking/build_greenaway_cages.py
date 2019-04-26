#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build cage structures in a brute force from precursors for greenaway paper.
Based off andrew_marsh_structures/build_cages.py.

Author: Andrew Tarzia

Date Created: 02 Apr 2019
"""

import glob
import stk
import sys
import os
import pandas as pd
sys.path.insert(0, '/home/atarzia/thesource/')
from stk_functions import topo_2_property, build_and_opt_cage, atarzia_MD_settings
from pywindow_functions import analyze_cage_from_MOL
from IO_tools import convert_PDB_2_XYZ


def check_done(NAME, output_csv):
    data = pd.read_csv(output_csv)
    if NAME not in list(set(data.name)):
        return False
    return True


def amine_defn(no):
    '''Output parameters for each amine.

    '''
    dict = {'A': (),
            'B': (),
            'C': (),
            }

    return dict[no]


def aldehyde_defn(no):
    '''Output parameters for each aldehyde.

    (list_no)
    '''
    dict = {'1': (1),
            '2': (1),
            '3': (1),
            '4': (1),
            '5': (1),
            '6': (1),
            '7': (1),
            '8': (1),
            '9': (1),
            '10': (1),
            '11': (2),
            '12': (2),
            '13': (2),
            '14': (2),
            '15': (2),
            '16': (2),
            '17': (2),
            '18': (2),
            '19': (2),
            '20': (2),
            '21': (2),
            '22': (3),
            '23': (3),
            '24': (3),
            '25': (3),
            '26': (3),
            }

    return dict[no]


def main():
    """Run script.

    """
    if (not len(sys.argv) == 1):
        print("""
Usage: build_cages.py output_file wipe run_build
    output_file: file to output results
    wipe: t/T if wipe output file
    run_build: t/T if you want to run the build stage""")
        sys.exit()
    else:
        pass

    macromod_ = '/home/atarzia/software/schrodinger_install'
    aldehyde_dir = '/home/atarzia/projects/method_benchmarking/greenaway_structures/aldehydes/'
    amine_dir = '/home/atarzia/projects/method_benchmarking/greenaway_structures/amines/'
    # get precursor files
    aldehyde_files = sorted(glob.glob(aldehyde_dir + '*.mol'))
    aldehyde_names = [i.replace(aldehyde_dir, '').rstrip('.mol') for i in aldehyde_files]
    amine_files = sorted(glob.glob(amine_dir + '*.mol'))
    amine_names = [i.replace(amine_dir, '').rstrip('.mol') for i in amine_files]

    # iterate over amines
    for i, amine in enumerate(amine_files):
        print(amine)
        bb_amine = stk.StructUnit3(amine, ['amine'])
        # iterate over aldehydes
        for j, aldehyde in enumerate(aldehyde_files):
            print(aldehyde)
            param = aldehyde_defn(aldehyde_names[j])
            lis = param  # [0]
            if lis == 1 or lis == 2:
                bb_alde = stk.StructUnit2(aldehyde, ['aldehyde'])
                topology_names = ['2p3', '4p6', '4p62', '6p9', 'dodec', '8p12']
            elif lis == 3:
                bb_alde = stk.StructUnit3(aldehyde, ['aldehyde'])
                topology_names = ['1p1', '4p4']  # , '2p2']
            topology_options = [topo_2_property(i, property='stk_func')
                                for i in topology_names]
            for k, topo in enumerate(topology_options):
                # naming convention: aldehyde-name_amine-name_topology
                NAME = aldehyde_names[j]
                NAME += '_' + amine_names[i] + '_'
                NAME += topology_names[k]
                prop_file = NAME + '_opt_properties.json'
                mole_file = NAME + '_opt_PWout.xyz'
                print('doing:', NAME)
                if os.path.isfile(NAME + '_opt.mol') is False:
                    # build cage and run optimization
                    cage = build_and_opt_cage(prefix=NAME,
                                              BB1=bb_alde,
                                              BB2=bb_amine,
                                              topology=topo,
                                              macromod_=macromod_,
                                              pdb=True,
                                              settings=atarzia_MD_settings())
                    # convert .pdb to .xyz using ASE
                    pdb = NAME + '_opt.pdb'
                    _, _ = convert_PDB_2_XYZ(pdb)
                    del _
                # check if completed and run pywindow if so
                if os.path.isfile(NAME + '_opt.mol') is True:
                    if os.path.isfile(prop_file) is False:
                        analyze_cage_from_MOL(file=NAME + '_opt.mol',
                                              prop_file=prop_file,
                                              mole_file=mole_file,
                                              include_coms=True)


if __name__ == "__main__":
    main()
