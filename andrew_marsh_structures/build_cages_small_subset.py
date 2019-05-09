#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build cage structures in a brute force way. Based off notebooks
build_cages_amines*.ipynb.

This script was used for the finer screening of top candidates. See ./build_cages.py
for the brute force, crude screening.

Author: Andrew Tarzia

Date Created: 30 Apr 2019
"""

import glob
import stk
import sys
import os
from rdkit.Chem import AllChem as Chem
sys.path.insert(0, '/home/atarzia/thesource/')
import andrew_marsh_structures.build_cages as BC
import IO_tools
import stk_functions
import pywindow_functions
import rdkit_functions


def precursor_pairings(amine_CN):
    '''Defined dictionary of precursor pairings

    '''
    if amine_CN == 2:
        # diamines
        pair_dict = {'aldehyde1': ['ami_5', 'ami_6', 'ami_7', 'ami_8eq',
                                   'ami_9eqtrans', 'ami_12', 'ami_13'],
                     'aldehyde2': ['ami_10eq', 'ami_11', 'ami_7', 'ami_14'],
                     'aldehyde3': ['ami_10eq', 'ami_11', 'ami_7', 'ami_14'],
                     'aldehyde4': ['ami_5', 'ami_6', 'ami_7', 'ami_8eq',
                                   'ami_9eqtrans', 'ami_12', 'ami_13']}

    elif amine_CN == 3:
        # triamines
        pair_dict = {'aldehyde1': ['ami_1ax', 'ami_3', 'ami_4'],
                     'aldehyde2': ['ami_1ax', 'ami_2ax'],
                     'aldehyde3': ['ami_1ax', 'ami_2ax'],
                     'aldehyde4': ['ami_1ax', 'ami_3', 'ami_4']}

    return pair_dict


def main():
    """Run script.

    """
    if (not len(sys.argv) == 3):
        print("""
Usage: build_cages.py output_file wipe run_build
    output_file: file to output results
    wipe: t/T if wipe output file""")
        sys.exit()
    else:
        output_file = sys.argv[1]
        wipe = sys.argv[2]

    macromod_ = '/home/atarzia/software/schrodinger_install'
    base_dir = '/home/atarzia/projects/andrew_marsh_structures/smaller_subset'
    alde_dir = os.path.join(base_dir, 'aldehydes/')
    ami2_dir = os.path.join(base_dir, 'diamines/')
    ami3_dir = os.path.join(base_dir, 'triamines/')
    # get precursor files
    alde_files = sorted(glob.glob(os.path.join(alde_dir, '*.mol')))
    print(alde_files)
    alde_names = [i.replace(alde_dir, '').rstrip('.mol') for i in alde_files]
    # read in mol files
    alde_struc = [stk.StructUnit3(i, ['aldehyde']) for i in alde_files]
    alde_smiles = [Chem.MolToSmiles(i.mol) for i in alde_struc]
    mol_list = [Chem.MolFromSmiles(i) for i in alde_smiles]
    mol_names = [i for i in alde_names]
    rdkit_functions.mol_list2grid(molecules=mol_list, names=mol_names,
                                  filename='aldehyde_precusors', mol_per_row=2,
                                  maxrows=10)

    # prepare output file
    output_csv = output_file
    if wipe.lower() == 't':
        if input('are you sure you wanna wipe? (t)') == 't':
            with open(output_csv, 'w') as f:
                f.write('name,bb1,SA1,bb2,SA2,topo,max_diam,p_diam,p_vol,p_diam_opt')
                f.write(',p_vol_opt,w_no,w_max,w_min,w_avg,w_diff,collapse,asym\n')

    # iterate over aldehydes
    for i, alde in enumerate(alde_files):
        # make diamines
        ami2_pairs = precursor_pairings(amine_CN=2)[alde_names[i]]
        print('aldehyde {} diamine pairs: {}'.format(alde, ami2_pairs))
        for ami2_name in ami2_pairs:
            ami2_file = os.path.join(ami2_dir, ami2_name+'.mol')
            ami2_struc = stk.StructUnit2(ami2_file, ['amine'])
            topology_names = ['2p3', '4p6', '4p62', '6p9']
            topology_options = [stk_functions.topo_2_property(i, property='stk_func')
                                for i in topology_names]
            for k, topo in enumerate(topology_options):
                # naming convention: aldehyde-name_amine-name_topology
                NAME = alde_names[i]
                NAME += '_' + ami2_name + '_'
                NAME += topology_names[k]
                prop_file = NAME + '_opt_properties.json'
                mole_file = NAME + '_opt_PWout.xyz'
                print('doing:', NAME)
                if os.path.isfile(NAME+'_opt.mol') is False:
                    # build cage and run optimization
                    cage = stk_functions.build_and_opt_cage(prefix=NAME,
                                                            BB1=alde_struc[i],
                                                            BB2=ami2_struc,
                                                            topology=topo,
                                                            macromod_=macromod_,
                                                            pdb=True,
                                                            settings=stk_functions.atarzia_long_MD_settings())
                    # convert .pdb to .xyz using ASE
                    pdb = NAME + '_opt.pdb'
                    _, _ = IO_tools.convert_PDB_2_XYZ(pdb)
                    del _
                # check if completed and run pywindow if so
                if os.path.isfile(NAME + '_opt.mol') is True:
                    if os.path.isfile(prop_file) is False:
                        pywindow_functions.analyze_cage_from_MOL(file=NAME+'_opt.mol',
                                                                 prop_file=prop_file,
                                                                 mole_file=mole_file,
                                                                 include_coms=True)
                # if pywindow is complete then analyse the cage and write out
                if os.path.isfile(prop_file) is True:
                    if BC.check_done(NAME, output_csv) is True:
                        continue
                    # reload cage if optimization is skipped
                    if 'cage' not in locals():
                        # cage does not exist
                        cage = stk.Cage([alde_struc[i], ami2_struc], topo)
                        cage.update_from_mol(NAME+'_opt.mol')
                    BC.assign_cage_properties(NAME=NAME, cage=cage,
                                              output_csv=output_csv)

        # make triamines
        ami3_pairs = precursor_pairings(amine_CN=3)[alde_names[i]]
        print('aldehyde {} triamine pairs: {}'.format(alde, ami3_pairs))
        for ami3_name in ami3_pairs:
            ami3_file = os.path.join(ami3_dir, ami3_name+'.mol')
            ami3_struc = stk.StructUnit3(ami3_file, ['amine'])
            topology_names = ['1p1', '4p4']
            topology_options = [stk_functions.topo_2_property(i, property='stk_func')
                                for i in topology_names]
            for k, topo in enumerate(topology_options):
                # naming convention: aldehyde-name_amine-name_topology
                NAME = alde_names[i]
                NAME += '_' + ami3_name + '_'
                NAME += topology_names[k]
                prop_file = NAME + '_opt_properties.json'
                mole_file = NAME + '_opt_PWout.xyz'
                print('doing:', NAME)
                if os.path.isfile(NAME+'_opt.mol') is False:
                    # build cage and run optimization
                    cage = stk_functions.build_and_opt_cage(prefix=NAME,
                                                            BB1=alde_struc[i],
                                                            BB2=ami3_struc,
                                                            topology=topo,
                                                            macromod_=macromod_,
                                                            pdb=True,
                                                            settings=stk_functions.atarzia_long_MD_settings())
                    # convert .pdb to .xyz using ASE
                    pdb = NAME + '_opt.pdb'
                    _, _ = IO_tools.convert_PDB_2_XYZ(pdb)
                    del _
                # check if completed and run pywindow if so
                if os.path.isfile(NAME + '_opt.mol') is True:
                    if os.path.isfile(prop_file) is False:
                        BC.analyze_cage_from_MOL(file=NAME+'_opt.mol',
                                                 prop_file=prop_file,
                                                 mole_file=mole_file,
                                                 include_coms=True)
                # if pywindow is complete then analyse the cage and write out
                if os.path.isfile(prop_file) is True:
                    if BC.check_done(NAME, output_csv) is True:
                        continue
                    # reload cage if optimization is skipped
                    if 'cage' not in locals():
                        # cage does not exist
                        cage = stk.Cage([alde_struc[i], ami3_struc], topo)
                        cage.update_from_mol(NAME+'_opt.mol')
                    BC.assign_cage_properties(NAME=NAME, cage=cage,
                                              output_csv=output_csv)
    # brute_analysis(output_csv, amine_type,
    #                precursor_dir, precursor_files, DB, amines,
    #                macromod_, des_topo, amine_SA)


if __name__ == "__main__":
    main()
