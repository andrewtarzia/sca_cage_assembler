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
import sys
import pandas as pd
import copy
import os
import matplotlib.pyplot as plt
sys.path.insert(0, '/home/atarzia/thesource/')
import andrew_marsh_structures.build_cages as BC
import stk_f


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


def subset_analysis(full_dataset, ami2_dir, ami3_dir):
    '''Function to run all analysis and screening.

    '''
    ald_colo = {'aldehyde1': 'k', 'aldehyde2': 'b',
                'aldehyde3': 'r', 'aldehyde4': 'g'}
    # full_dataset = pd.read_csv(output_csv)
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
    if (len(sys.argv) < 2):
        logging.info(f"""
Usage: smaller_subset_analysis.py *output_files
    output_file (str): CSV files with build results. Merges all CSVs listed""")
        sys.exit()
    else:
        output_files = sys.argv[1:]

    print(output_files)
    DBs = [pd.read_csv(i) for i in output_files[1:]]
    print(DBs)
    full_dataset = pd.read_csv(output_files[0])
    for i in DBs:
        full_dataset = full_dataset.append(i)

    print(full_dataset)
    # sys.exit()

    base_dir = '/home/atarzia/projects/andrew_marsh_structures/smaller_subset'
    ami2_dir = os.path.join(base_dir, 'diamines/')
    ami3_dir = os.path.join(base_dir, 'triamines/')

    # run analysis for both amine types
    subset_analysis(full_dataset=full_dataset,
                    ami2_dir=ami2_dir, ami3_dir=ami3_dir)


if __name__ == "__main__":
    # logging.basicConfig(level=logging.DEBUG, format='%(levelname)s-%(message)s')
    # logging.debug(f'Debug mode!')
    logging.basicConfig(level=logging.INFO, format='')
    main()
