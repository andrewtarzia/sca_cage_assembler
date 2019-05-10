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
import pandas as pd
import copy
import os
import matplotlib.pyplot as plt
from rdkit.Chem import AllChem as Chem
sys.path.insert(0, '/home/atarzia/thesource/')
import andrew_marsh_structures.build_cages as BC
import IO_tools
import stk_f
import pywindow_f
import rdkit_f


def plots_amines2(final_db, ald_colo, DB, des_topo):
    #########################################################################
    fig, ax = plt.subplots()
    for i in ald_colo:
        data = final_db[final_db.bb1 == i]
        ax.scatter(data.w_max, data.p_diam_opt, c=ald_colo[i], alpha=0.8,
                   edgecolor='k', marker='o', s=80, label=i)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('max window size [$\mathrm{\AA}$]', fontsize=16)
    ax.set_ylabel('pore diameter [$\mathrm{\AA}$]', fontsize=16)
    ax.set_xlim(0, 40)
    ax.set_ylim(0, 40)
    ax.legend(fontsize=12)
    fig.tight_layout()
    fig.savefig("pdiam_VS_wmax_2f.pdf",
                dpi=720, bbox_inches='tight')
    plt.close()
    #########################################################################
    for topo in des_topo:
        temp_db = final_db[final_db.topo == topo]
        fig, ax = plt.subplots()
        for i in ald_colo:
            data = temp_db[temp_db.bb1 == i]
            ax.scatter(data.w_max, data.p_diam_opt, c=ald_colo[i], alpha=0.8,
                       edgecolor='k', marker='o', s=80, label=i)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('max window size [$\mathrm{\AA}$]', fontsize=16)
        ax.set_ylabel('pore diameter [$\mathrm{\AA}$]', fontsize=16)
        ax.set_xlim(0, 40)
        ax.set_ylim(0, 40)
        ax.legend(fontsize=12)
        fig.tight_layout()
        fig.savefig("pdiam_VS_wmax_2f_"+topo+".pdf",
                    dpi=720, bbox_inches='tight')
        plt.close()
    #########################################################################
    for i in ald_colo:
        fig, ax = plt.subplots()
        data = final_db[final_db.bb1 == i]
        ax.scatter(data.asym, data.p_diam_opt, c=ald_colo[i],
                   alpha=0.8,
                   edgecolor='k', marker='o', s=80, label=i)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('average asymmetry in window size [$\mathrm{\AA}$]',
                      fontsize=16)
        ax.set_ylabel('pore diameter [$\mathrm{\AA}$]', fontsize=16)
        ax.set_xlim(0, 5)
        ax.set_ylim(0, 40)
        ax.legend(fontsize=12)
        fig.tight_layout()
        fig.savefig("pdiam_VS_asym_2f_"+i+".pdf",
                    dpi=720, bbox_inches='tight')
        plt.close()
    #########################################################################
    for a in ald_colo:
        data = final_db[final_db.bb1 == a]
        filename = 'candidate_amines2f_'+a+'_largepore'
        BC.output_precursor_struct(data, filename, DB,
                                   sorter='p_diam_opt', rev=True,
                                   prop1='p_diam_opt', prop2='p_diam_opt')
    #########################################################################


def plots_amines3(final_db, ald_colo, DB):
    #########################################################################
    fig, ax = plt.subplots()
    for i in ald_colo:
        data = final_db[final_db.bb1 == i]
        ax.scatter(data.w_max, data.p_diam_opt, c=ald_colo[i], alpha=0.8,
                   edgecolor='k', marker='o', s=80, label=i)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('max window size [$\mathrm{\AA}$]', fontsize=16)
    ax.set_ylabel('pore diameter [$\mathrm{\AA}$]', fontsize=16)
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 20)
    ax.legend(fontsize=12)
    fig.tight_layout()
    fig.savefig("pdiam_VS_wmax_3f.pdf",
                dpi=720, bbox_inches='tight')
    plt.close()
    #########################################################################
    temp_db = final_db[final_db.topo == '1p1']
    fig, ax = plt.subplots()
    for i in ald_colo:
        data = temp_db[temp_db.bb1 == i]
        ax.scatter(data.w_max, data.p_diam_opt, c=ald_colo[i], alpha=0.8,
                   edgecolor='k', marker='o', s=80, label=i)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('max window size [$\mathrm{\AA}$]', fontsize=16)
    ax.set_ylabel('pore diameter [$\mathrm{\AA}$]', fontsize=16)
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 20)
    ax.legend(fontsize=12)
    fig.tight_layout()
    fig.savefig("pdiam_VS_wmax_3f_1p1.pdf",
                dpi=720, bbox_inches='tight')
    plt.close()
    #########################################################################
    temp_db = final_db[final_db.topo == '4p4']
    fig, ax = plt.subplots()
    for i in ald_colo:
        data = temp_db[temp_db.bb1 == i]
        ax.scatter(data.w_max, data.p_diam_opt, c=ald_colo[i], alpha=0.8,
                   edgecolor='k', marker='o', s=80, label=i)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('max window size [$\mathrm{\AA}$]', fontsize=16)
    ax.set_ylabel('pore diameter [$\mathrm{\AA}$]', fontsize=16)
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 20)

    ax.legend(fontsize=12)
    fig.tight_layout()
    fig.savefig("pdiam_VS_wmax_3f_4p4.pdf",
                dpi=720, bbox_inches='tight')
    plt.close()
    #########################################################################
    for i in ald_colo:
        fig, ax = plt.subplots()
        data = final_db[final_db.bb1 == i]
        ax.scatter(data.w_max, data.p_diam_opt, c=ald_colo[i], alpha=0.8,
                   edgecolor='k', marker='o', s=80, label=i)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('max window size [$\mathrm{\AA}$]', fontsize=16)
        ax.set_ylabel('pore diameter [$\mathrm{\AA}$]', fontsize=16)
        ax.set_xlim(0, 20)
        ax.set_ylim(0, 20)
        ax.legend(fontsize=12)
        fig.tight_layout()
        fig.savefig("pdiam_VS_wmax_3f_"+i+".pdf",
                    dpi=720, bbox_inches='tight')
        plt.close()
    #########################################################################
    final_2_db = final_db[final_db.topo == '4p4']
    for i in ald_colo:
        fig, ax = plt.subplots()
        data = final_2_db[final_2_db.bb1 == i]
        ax.scatter(data.asym, data.p_diam_opt, c=ald_colo[i],
                   alpha=0.8,
                   edgecolor='k', marker='o', s=80, label=i)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('average asymmetry in window size [$\mathrm{\AA}$]',
                      fontsize=16)
        ax.set_ylabel('pore diameter [$\mathrm{\AA}$]', fontsize=16)
        ax.set_xlim(0, 5)
        ax.set_ylim(0, 20)
        ax.legend(fontsize=12)
        fig.tight_layout()
        fig.savefig("pdiam_VS_asym_3f_"+i+".pdf",
                    dpi=720, bbox_inches='tight')
        plt.close()
    #########################################################################
    for a in ald_colo:
        data = final_db[final_db.bb1 == a]
        filename = 'candidate_amines3f_'+a
        BC.output_precursor_struct(data, filename, DB,
                                   sorter='p_diam_opt',
                                   prop1='p_diam_opt', prop2='p_diam_opt')
    #########################################################################


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


def screening_process(dataset):
    '''Screen cages based on some filters in this function.

    Returns new DataFrame

    '''
    final_db = pd.DataFrame(columns=dataset.columns)
    # iterate through and collect amine buidling blocks that produce cages
    # with reasonable properties
    for i, row in dataset.iterrows():
        ################################
        # remove certain topologies
        if row.topo in ['4p62']:
            continue
        ################################
        # remove those with pore diamter < 3.4 angstrom
        # (Computationally-inspired discovery of an unsymmetrical porous
        # organic cage)
        if row.p_diam_opt < 3.4:
            continue
        ################################
        # remove those that have no windows
        if row.w_no == 0:
            continue
        ################################
        # remove structures with no. windows < expected for topology
        if row.w_no < stk_f.topo_2_property(row.topo, property='expected_wind'):
            continue
        ################################
        # remove those with max window diamter < 2.8 angstrom
        # (Computationally-inspired discovery of an unsymmetrical porous
        # organic cage)
        if row.w_max < 2.8:
            continue
        ################################
        # recalculate asymetry - remove cases with asymetry > XX
        # do not use the asymetry defined ni stk_f (deprecated)
        # It does not handle different window types. Use stk window_difference
        # function
        asymmetry = row.w_diff
        # no windows found
        if asymmetry == 'None':
            continue
        row.asym = float(asymmetry)
        if row.asym < 0:
            continue
        if row.asym > 2:
            continue
        final_db = final_db.append(row)
    return final_db


def subset_analysis(output_csv, ami2_dir, ami3_dir):
    '''Function to run all analysis and screening.

    '''
    ald_colo = {'aldehyde1': 'k', 'aldehyde2': 'b',
                'aldehyde3': 'r', 'aldehyde4': 'g'}
    full_dataset = pd.read_csv(output_csv)
    working_dataset = copy.deepcopy(full_dataset)
    # do not get cage energies here
    # screen cages
    final_db = screening_process(dataset=working_dataset)
    logging.info(f'>>> {len(final_db)} cages remaining after screening!')
    logging.info(f'>>> doing all plotting')
    # split into di and tri amines
    diamine_db = pd.DataFrame(columns=final_db.columns)
    triamine_db = pd.DataFrame(columns=final_db.columns)
    logging.info(f'{len(diamine_db)} + {len(triamine_db)} = {len(final_db)}')
    for i, row in final_db.iterrows():
        if row.topo in ['2p3', '4p6', '6p9']:
            diamine_db = diamine_db.append(row)
        elif row.topo in ['1p1', '4p4']:
            triamine_db = triamine_db.append(row)
    logging.info(f'{len(diamine_db)} + {len(triamine_db)} = {len(final_db)}')
    # diamines
    des_topo = ['2p3', '4p6', '6p9']
    plots_amines2(final_db=diamine_db, ald_colo=ald_colo, DB=ami2_dir,
                  des_topo=des_topo)
    # triamines
    plots_amines3(final_db=triamine_db, ald_colo=ald_colo, DB=ami3_dir)
    logging.info(f'>>> done!')


def main():
    """Run script.

    """
    if (not len(sys.argv) == 3):
        logging.info(f"""
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
    if wipe.lower() == 't':
        if input('are you sure you wanna wipe? (t)') == 't':
            with open(output_csv, 'w') as f:
                f.write('name,bb1,SA1,bb2,SA2,topo,max_diam,p_diam,p_vol,p_diam_opt')
                f.write(',p_vol_opt,w_no,w_max,w_min,w_avg,w_diff,collapse,asym\n')

    # iterate over aldehydes
    for i, alde in enumerate(alde_files):
        # make diamines
        ami2_pairs = precursor_pairings(amine_CN=2)[alde_names[i]]
        logging.info(f'aldehyde {alde} diamine pairs: {ami2_pairs}')
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
                    cage = stk_f.build_and_opt_cage(prefix=NAME,
                                                            BB1=alde_struc[i],
                                                            BB2=ami2_struc,
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
                        cage = stk.Cage([alde_struc[i], ami2_struc], topo)
                        cage.update_from_mol(NAME+'_opt.mol')
                    BC.assign_cage_properties(NAME=NAME, cage=cage,
                                              output_csv=output_csv)

        # make triamines
        ami3_pairs = precursor_pairings(amine_CN=3)[alde_names[i]]
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
    # run analysis for both amine types
    subset_analysis(output_csv, ami2_dir=ami2_dir, ami3_dir=ami3_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, format='%(levelname)s-%(message)s')
    logging.debug(f'Debug mode!')
    # logging.basicConfig(level=logging.INFO, format='')
    main()
