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

import logging
import glob
import stk
import sys
import os
from rdkit.Chem import AllChem as Chem
sys.path.insert(0, '/home/atarzia/thesource/')
import andrew_marsh_structures.build_cages as BC
import IO_tools
import stk_f
import pywindow_f
import rdkit_f


def precursor_pairings(amine_CN, aldehyde):
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
    if aldehyde in pair_dict:
        return pair_dict[aldehyde]
    else:
        return []


def main():
    """Run script.

    """
    if (not len(sys.argv) == 4):
        logging.info(f"""
Usage: build_cages_small_subset.py output_file wipe targ_aldehyde
    output_file (str): file to output results
    wipe (str): t/T if wipe output file
    targ_aldehyde (str): name of aldehyde to build/run analysis on""")
        sys.exit()
    else:
        output_file = sys.argv[1]
        wipe = sys.argv[2]
        targ_aldehyde = sys.argv[3]

    macromod_ = '/home/atarzia/software/schrodinger_install'
    base_dir = '/home/atarzia/projects/andrew_marsh_structures/smaller_subset'
    macromod_output = os.path.join(os.getcwd(), 'data/')
    alde_dir = os.path.join(base_dir, 'aldehydes/')
    ami2_dir = os.path.join(base_dir, 'diamines/')
    ami3_dir = os.path.join(base_dir, 'triamines/')
    # get precursor files
    alde_files = [i for i in sorted(glob.glob(os.path.join(alde_dir, '*.mol')))
                  if '_opt.mol' not in i]
    logging.debug(alde_files)
    alde_names = [i.replace(alde_dir, '').rstrip('.mol') for i in alde_files]
    # read in mol files
    alde_struc = [stk.StructUnit3(i, ['aldehyde']) for i in alde_files]
    alde_smiles = [Chem.MolToSmiles(i.mol) for i in alde_struc]
    mol_list = [Chem.MolFromSmiles(i) for i in alde_smiles]
    mol_names = [i for i in alde_names]
    rdkit_f.mol_list2grid(molecules=mol_list, names=mol_names,
                                  filename='aldehyde_precusors', mol_per_row=2,
                                  maxrows=10)

    # prepare output file
    output_csv = output_file
    if os.path.isfile(output_csv) is False:
        with open(output_csv, 'w') as f:
            f.write('name,bb1,SA1,bb2,SA2,topo,max_diam,p_diam,p_vol,p_diam_opt')
            f.write(',p_vol_opt,w_no,w_max,w_min,w_avg,w_diff,collapse,asym\n')
    if wipe.lower() == 't':
        if input('are you sure you wanna wipe? (t)') == 't':
            with open(output_csv, 'w') as f:
                f.write('name,bb1,SA1,bb2,SA2,topo,max_diam,p_diam,p_vol,p_diam_opt')
                f.write(',p_vol_opt,w_no,w_max,w_min,w_avg,w_diff,collapse,asym\n')

    # iterate over aldehydes
    for i, alde in enumerate(alde_files):
        if alde_names[i] != targ_aldehyde:
            continue
        # make diamines
        ami2_pairs = precursor_pairings(amine_CN=2, aldehyde=alde_names[i])
        logging.info(f'aldehyde {alde_names[i]} diamine pairs: {ami2_pairs}')
        for ami2_name in ami2_pairs:
            ami2_file = os.path.join(ami2_dir, ami2_name+'.mol')
            ami2_struc = stk.StructUnit2(ami2_file, ['amine'])
            topology_names = ['2p3', '4p6', '4p62', '6p9']
            topology_options = [stk_f.topo_2_property(i, property='stk_func')
                                for i in topology_names]
            for k, topo in enumerate(topology_options):
                # naming convention: aldehyde-name_amine-name_topology
                NAME = alde_names[i]
                NAME += '_' + ami2_name + '_'
                NAME += topology_names[k]
                prop_file = NAME + '_opt_properties.json'
                mole_file = NAME + '_opt_PWout.xyz'
                logging.info(f'doing:{NAME}')
                if os.path.isfile(NAME+'_opt.mol') is False:
                    # build cage and run optimization
                    sim_output = macromod_output + NAME
                    cage = stk_f.build_and_opt_cage(prefix=NAME,
                                                    BB1=alde_struc[i],
                                                    BB2=ami2_struc,
                                                    topology=topo,
                                                    macromod_=macromod_,
                                                    pdb=True,
                                                    settings=stk_f.atarzia_long_MD_settings(),
                                                    output_dir=sim_output)
                    # convert .pdb to .xyz using ASE
                    pdb = NAME + '_opt.pdb'
                    _, _ = IO_tools.convert_PDB_2_XYZ(pdb)
                    del _
                # check if completed and run pywindow if so
                if os.path.isfile(NAME + '_opt.mol') is True:
                    if os.path.isfile(prop_file) is False:
                        pywindow_f.analyze_cage_from_MOL(file=NAME+'_opt.mol',
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
        ami3_pairs = precursor_pairings(amine_CN=3, aldehyde=alde_names[i])
        logging.info(f'aldehyde {alde} triamine pairs: {ami3_pairs}')
        for ami3_name in ami3_pairs:
            ami3_file = os.path.join(ami3_dir, ami3_name+'.mol')
            ami3_struc = stk.StructUnit3(ami3_file, ['amine'])
            topology_names = ['1p1', '4p4']
            topology_options = [stk_f.topo_2_property(i, property='stk_func')
                                for i in topology_names]
            for k, topo in enumerate(topology_options):
                # naming convention: aldehyde-name_amine-name_topology
                NAME = alde_names[i]
                NAME += '_' + ami3_name + '_'
                NAME += topology_names[k]
                prop_file = NAME + '_opt_properties.json'
                mole_file = NAME + '_opt_PWout.xyz'
                logging.info(f'doing:{NAME}')
                if os.path.isfile(NAME+'_opt.mol') is False:
                    # build cage and run optimization
                    cage = stk_f.build_and_opt_cage(prefix=NAME,
                                                            BB1=alde_struc[i],
                                                            BB2=ami3_struc,
                                                            topology=topo,
                                                            macromod_=macromod_,
                                                            pdb=True,
                                                            settings=stk_f.atarzia_long_MD_settings())
                    # convert .pdb to .xyz using ASE
                    pdb = NAME + '_opt.pdb'
                    _, _ = IO_tools.convert_PDB_2_XYZ(pdb)
                    del _
                # check if completed and run pywindow if so
                if os.path.isfile(NAME + '_opt.mol') is True:
                    if os.path.isfile(prop_file) is False:
                        pywindow_f.analyze_cage_from_MOL(file=NAME+'_opt.mol',
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


if __name__ == "__main__":
    # logging.basicConfig(level=logging.DEBUG, format='%(levelname)s-%(message)s')
    # logging.debug(f'Debug mode!')
    logging.basicConfig(level=logging.INFO, format='')
    main()
