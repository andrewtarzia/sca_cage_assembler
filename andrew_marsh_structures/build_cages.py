#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build cage structures in a brute force way. Based off notebooks
build_cages_amines*.ipynb.

Author: Andrew Tarzia

Date Created: 18 Mar 2019
"""

import glob
import stk
import sys
import os
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import copy
from stk.molecular.molecules import MacroMoleculeBuildError
sys.path.insert(0, '/home/atarzia/thesource/')
import stk_functions
import pywindow_functions


def output_precursor_struct(data, filename, DB, prop1, prop2, sorter,
                            rev=False):
    '''Produce figure ranking precursors by 'sorter'

    '''
    smiles_done = []
    max_c = 9
    if rev is False:
        sorted_data = data.sort_values(by=[sorter])
    else:
        sorted_data = data.sort_values(by=[sorter], ascending=False)
    sorted_mols = []
    sorted_props = []
    sorted_NAMES = []
    for i, row in sorted_data.iterrows():
        NAME = DB+row.bb2+'.mol'
        smiles = Chem.MolToSmiles(Chem.MolFromMolFile(NAME))
        if smiles not in smiles_done:
            sorted_NAMES.append(row.bb1+'_'+row.bb2+'_'+row.topo)
            smiles_done.append(smiles)
            sorted_mols.append(Chem.MolFromSmiles(smiles))
            sorted_props.append((row[prop1], row[prop2]))
    ranges = np.arange(0, len(sorted_mols)+1, max_c)
    # print(ranges)
    for i in ranges:
        if i > max_c:
            break
        # print(a, i)
        _set = sorted_mols[i: i+max_c]
        _set_props = sorted_props[i: i+max_c]
        if len(_set) == 0:
            break
        print(sorted_NAMES[i: i+max_c])
        img = Draw.MolsToGridImage(_set, molsPerRow=3,
                                   subImgSize=(125, 125),
                                   # legends=[str(round(i[0], 3))+' ('+str(round(i[1], 2))+')'
                                   legends=[str(round(i[1], 2))
                                            for i in _set_props],
                                   useSVG=False)
        img.save(filename+'_'+str(i)+'.png')


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
        ax.scatter(data.SA2, data.p_diam_opt, c=ald_colo[i], alpha=0.8,
                   edgecolor='k', marker='o', s=80, label=i)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('SAscore of amine', fontsize=16)
        ax.set_ylabel('pore diameter [$\mathrm{\AA}$]', fontsize=16)
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 40)
        ax.legend(fontsize=12)
        fig.tight_layout()
        fig.savefig("pdiam_VS_SA2_2f_"+i+".pdf",
                    dpi=720, bbox_inches='tight')
        plt.close()
    #########################################################################
    for i in ald_colo:
        fig, ax = plt.subplots()
        data = final_db[final_db.bb1 == i]
        ax.scatter(data.assym/data.w_no, data.p_diam_opt, c=ald_colo[i],
                   alpha=0.8,
                   edgecolor='k', marker='o', s=80, label=i)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('average asymmetry in window size [$\mathrm{\AA}$]', fontsize=16)
        ax.set_ylabel('pore diameter [$\mathrm{\AA}$]', fontsize=16)
        ax.set_xlim(0, 20)
        ax.set_ylim(0, 40)
        ax.legend(fontsize=12)
        fig.tight_layout()
        fig.savefig("pdiam_VS_asym_2f_"+i+".pdf",
                    dpi=720, bbox_inches='tight')
        plt.close()
    #########################################################################
    for a in ald_colo:
        data = final_db[final_db.bb1 == a]
        filename = 'candidate_amines2f_'+a
        output_precursor_struct(data, filename, DB,
                                sorter='SA2',
                                prop1='SA2', prop2='p_diam_opt')
    #########################################################################
    for a in ald_colo:
        data = final_db[final_db.bb1 == a]
        filename = 'candidate_amines2f_'+a+'_largepore'
        output_precursor_struct(data, filename, DB,
                                sorter='p_diam_opt', rev=True,
                                prop1='SA2', prop2='p_diam_opt')
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
        ax.scatter(data.SA2, data.p_diam_opt, c=ald_colo[i], alpha=0.8,
                   edgecolor='k', marker='o', s=80, label=i)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('SAscore of amine', fontsize=16)
        ax.set_ylabel('pore diameter [$\mathrm{\AA}$]', fontsize=16)
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 20)
        ax.legend(fontsize=12)
        fig.tight_layout()
        fig.savefig("pdiam_VS_SA2_3f_"+i+".pdf",
                    dpi=720, bbox_inches='tight')
        plt.close()
    #########################################################################
    for i in ald_colo:
        fig, ax = plt.subplots()
        data = final_2_db[final_2_db.bb1 == i]
        ax.scatter(data.assym/data.w_no, data.p_diam_opt, c=ald_colo[i],
                   alpha=0.8,
                   edgecolor='k', marker='o', s=80, label=i)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('average asymmetry in window size [$\mathrm{\AA}$]', fontsize=16)
        ax.set_ylabel('pore diameter [$\mathrm{\AA}$]', fontsize=16)
        ax.set_xlim(0, 12)
        ax.set_ylim(0, 20)
        ax.legend(fontsize=12)
        fig.tight_layout()
        fig.savefig("pdiam_VS_asym_3f_"+i+".pdf",
                    dpi=720, bbox_inches='tight')
        plt.close()
    #########################################################################
    fig, ax = plt.subplots()
    for i in ald_colo:
        data = final_2_db[final_2_db.bb1 == i]
        ax.scatter(data.cage_ey, data.p_diam_opt, c=ald_colo[i], alpha=0.8,
                   edgecolor='k', marker='o', s=80, label=i)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('cage energy [kJ/mol]', fontsize=16)
    ax.set_ylabel('pore diameter [$\mathrm{\AA}$]', fontsize=16)
    ax.set_xlim(-500, 5000)
    ax.set_ylim(0, 20)

    ax.legend(fontsize=12)
    fig.tight_layout()
    fig.savefig("pdiam_VS_cageE_3f.pdf",
                dpi=720, bbox_inches='tight')
    plt.close()
    #########################################################################
    fig, ax = plt.subplots()
    for i in ald_colo:
        data = final_2_db[final_2_db.bb1 == i]
        ax.scatter(data.form_ey, data.p_diam_opt, c=ald_colo[i], alpha=0.8,
                   edgecolor='k', marker='o', s=80, label=i)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('formation energy [kJ/mol]', fontsize=16)
    ax.set_ylabel('pore diameter [$\mathrm{\AA}$]', fontsize=16)
    ax.set_xlim(-3090, 3000)
    ax.set_ylim(0, 20)
    ax.legend(fontsize=12)
    fig.tight_layout()
    fig.savefig("pdiam_VS_formE_3f.pdf",
                dpi=720, bbox_inches='tight')
    plt.close()
    #########################################################################
    for i in ald_colo:
        fig, ax = plt.subplots()
        data = final_2_db[final_2_db.bb1 == i]
        ax.scatter(data.form_ey, data.p_diam_opt, c=ald_colo[i], alpha=0.8,
                   edgecolor='k', marker='o', s=80, label=i)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('formation energy [kJ/mol]', fontsize=16)
        ax.set_ylabel('pore diameter [$\mathrm{\AA}$]', fontsize=16)
        ax.set_xlim(-3090, 3000)
        ax.set_ylim(0, 20)
        ax.legend(fontsize=12)
        fig.tight_layout()
        fig.savefig("pdiam_VS_formE_3f_"+i+".pdf",
                    dpi=720, bbox_inches='tight')
        plt.close()
    #########################################################################
    for a in ald_colo:
        data = final_db[final_db.bb1 == a]
        filename = 'candidate_amines3f_'+a
        output_precursor_struct(data, filename, DB,
                                sorter='SA2',
                                prop1='SA2', prop2='p_diam_opt')
    #########################################################################


def check_output_numbers(precursor_dir, DB, topology_options, prefix):
    # Check that the number output files matches the expected number based on
    # inputs
    print(len(glob.glob(precursor_dir+'*mol')), 'precusrors in', precursor_dir)
    print(len(glob.glob(DB+'*mol')), 'precusrors in', DB)
    print(len(topology_options), 'possible topologies per cage')
    n_cages = len(glob.glob(precursor_dir+'*mol'))
    n_cages *= len(glob.glob(DB+'*mol')) * len(topology_options)
    print('>>', n_cages, 'cages')
    print(len(glob.glob('*'+prefix+'*_opt*mol')), 'cages optimized')
    # Check that the number output files matches expected number based on
    # inputs
    opt_cages = len(glob.glob('*'+prefix+'*_opt*mol'))
    print('>>', opt_cages, 'cages')
    print(len(glob.glob('*_properties.json')),
          'cage properties calculated')
    print(len(glob.glob('*_PWout.xyz')),
          'cage window structures output')


def check_done(NAME, output_csv):
    data = pd.read_csv(output_csv)
    if NAME not in list(set(data.name)):
        return False
    return True


def brute_cage_build(precursor_struc, precursor_names, precursor_files,
                     topology_names, topology_options, amines, DB, amine_type,
                     output_csv, macromod_):
    '''Nested loops for building all possible cages.

    '''
    no_amines = len(amines)
    for i, prec in enumerate(precursor_struc):
        for j, amine in enumerate(amines):
            if amine_type == '2':
                bb_amine = stk.StructUnit2(amine)
            elif amine_type == '3':
                bb_amine = stk.StructUnit3(amine)
            for k, topo in enumerate(topology_options):
                # naming convention: aldehyde-name_amine-name_topology
                # amine-name: amineDB_NO
                NAME = precursor_names[i].replace('.mol', '')
                NAME += '_'+amine.replace('.mol', '').replace(DB, '')+'_'
                NAME += topology_names[k]
                prop_file = NAME+'_opt_properties.json'
                mole_file = NAME+'_opt_PWout.xyz'
                print('doing:', NAME, ' -- amine:', j, 'of', no_amines)
                if os.path.isfile(NAME+'_opt.mol') is False:
                    # run optimization
                    try:
                        cage = stk.Cage([prec, bb_amine], topo)
                        cage.write(NAME+'.mol')
                        stk.macromodel_opt(cage, macromodel_path=macromod_,
                                           settings={'restricted': False,
                                                     'timeout': None,
                                                     'force_field': 16,
                                                     'max_iter': 2500,
                                                     'gradient': 0.05,
                                                     'md': True
                                                    },
                                           md={'timeout': None,
                                               'force_field': 16,
                                               'temp': 700,
                                               'confs': 50,
                                               'time_step': 1.0,
                                               'eq_time': 10,
                                               'sim_time': 200,
                                               'max_iter': 2500,
                                               'gradient': 0.05
                                              })
                        cage.write(NAME+'_opt.mol')
                    except MacroMoleculeBuildError:
                        pass
                # check if completed and run pywindow if so
                if os.path.isfile(NAME+'_opt.mol') is True:
                    if os.path.isfile(prop_file) is False:
                        pywindow_functions.run_on_cage(file=NAME+'_opt.mol',
                                                       prop_file=prop_file,
                                                       mole_file=mole_file)
                # if pywindow is complete then analyse the cage and write out
                if os.path.isfile(prop_file) is True:
                    if check_done(NAME, output_csv) is True:
                        continue
                    # reload cage if optimization is skipped
                    if 'cage' not in locals():
                        # cage does not exist
                        cage = stk.Cage([prec, bb_amine], topo)
                        cage.update_from_mol(NAME+'_opt.mol')
                    # all analysis done successfully -- output
                    bb1, SA1, bb2, SA2, topo = '-', '-', '-', '-', '-'
                    max_diam, p_diam, p_vol = '-', '-', '-'
                    p_diam_opt, p_vol_opt = '-', '-'
                    w_no, w_max, w_min, w_avg = '-', '-', '-', '-'
                    w_diff, collapsed, assymetry = '-', '-', '-'
                    with open(prop_file, 'r') as f:
                        data = json.load(f)
                    bb1 = NAME.split('_')[0]
                    bb2 = '_'.join(NAME.split('_')[1:3])
                    topo = NAME.split('_')[3]
                    max_diam = str(data['maximum_diameter']['diameter'])
                    p_diam = str(data['pore_diameter']['diameter'])
                    p_diam_opt = str(data['pore_diameter_opt']['diameter'])
                    p_vol = str(data['pore_volume'])
                    p_vol_opt = str(data['pore_volume_opt'])
                    if data['windows']['diameters'] is None:
                        w_no, w_max, w_min, w_avg, w_diff = '0', '0', '0', '0', '0'
                        collapse = '2'  # unsure
                        assymetry = '-2'
                    elif len(data['windows']['diameters']) == 0:
                        w_no, w_max, w_min, w_avg, w_diff = '0', '0', '0', '0', '0'
                        collapse = '1'  # collapsed
                        assymetry = '-2'
                    else:
                        if max(data['windows']['diameters']) < 200:
                            w_no = str(len(data['windows']['diameters']))
                            w_max = str(max(data['windows']['diameters']))
                            w_min = str(min(data['windows']['diameters']))
                            w_avg = str(np.average(data['windows']['diameters']))
                            if topo != '4p62':
                                w_diff = str(cage.window_difference())
                            else:
                                w_diff = None
                            if w_diff is None:
                                w_diff = '-5'
                            assymetry = str(stk_functions.get_asymmetry(data))
                            coll_flag = stk_functions.is_collapse(
                                topo=topo, avg_diff=w_diff,
                                max_window_diam=w_max,
                                cavity_size=cage.cavity_size(),
                                no_window=w_no)
                            if coll_flag is True:
                                collapse = '0'
                            elif coll_flag is False:
                                collapse = '1'
                            elif coll_flag is None:
                                collapse = '2'
                        else:
                            w_no, w_max, w_min, w_avg, w_diff = '-1', '-1', '-1', '-1', '-1'
                            collapse = '2'  # unsure
                            assymetry = '-1'
                    with open(output_csv, 'a') as f:
                        f.write(NAME+','+bb1+','+SA1+','+bb2+','+SA2+','+topo+',')
                        f.write(max_diam+','+p_diam+','+p_vol+','+p_diam_opt+',')
                        f.write(p_vol_opt+','+w_no+','+w_max+','+w_min+','+w_avg+',')
                        f.write(w_diff+','+collapse+','+assymetry)
                        f.write('\n')


def get_OPLS3_energy_of_list(out_file, structures, macromod_, dir='', opt=False):
    '''Get OPLS3 single point energy of a list of structures.

    Keyword Arguments:
        file (str) - file that tracks the structures done to avoid redoing
        structures (list) - list of structure files
        dir (str) - directory where structures are if not in working dir
        opt (bool) - True if optimization is required using macromodel
    '''
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
                               settings={'restricted': False,
                                         'timeout': None,
                                         'force_field': 16,
                                         'max_iter': 2500,
                                         'gradient': 0.05,
                                         'md': True},
                               md={'timeout': None,
                                   'force_field': 16,
                                   'temp': 700,
                                   'confs': 50,
                                   'time_step': 1.0,
                                   'eq_time': 10,
                                   'sim_time': 200,
                                   'max_iter': 2500,
                                   'gradient': 0.05})
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


def get_formation_energies(data):
    '''Calculate formation energies based on prod - react energies.

    '''
    # by definition
    H2O_energy = 0.0
    form_energies = []
    stoich = {'1p1': {'bb1': 1, 'bb2': 1, 'h2o': 3},
              '4p4': {'bb1': 4, 'bb2': 4, 'h2o': 12},
              '2p3': {'bb1': 2, 'bb2': 3, 'h2o': 6},
              '4p6': {'bb1': 4, 'bb2': 6, 'h2o': 12},
              '4p62': {'bb1': 4, 'bb2': 6, 'h2o': 12},
              '6p9': {'bb1': 6, 'bb2': 9, 'h2o': 18}}
    for i, row in data.iterrows():
        FE = (row.cage_ey + stoich[row.topo]['h2o'] * H2O_energy) - (row.bb1_ey * stoich[row.topo]['bb1'] + row.bb2_ey * stoich[row.topo]['bb2'])
        form_energies.append(FE)
    return form_energies


def screening_cages(dataset, des_topo, SA_data):
    '''Screen cages based on some filters in this function.

    Returns new DataFrame

    '''
    amines = []
    cages = []
    final_db = pd.DataFrame(columns=dataset.columns)
    # iterate through and collect amine buidling blocks that produce cages
    # with reasonable properties
    for i, row in dataset.iterrows():
        if row.topo not in des_topo:
            continue
        # remove those with pore diamter < 3.4 angstrom
        # (Computationally-inspired discovery of an unsymmetrical porous
        # organic cage)
        if row.p_diam_opt < 3.4:
            continue
        # remove those that have no windows
        if row.assym == -2:
            continue
        # remove those with max window diamter < 2.8 angstrom
        # (Computationally-inspired discovery of an unsymmetrical porous
        # organic cage)
        if row.p_diam_opt < 2.8:
            continue
        # recalculate assymetry - remove cases with assymetry > 1
        NAME = row.bb1+'_'+row.bb2+'_'+row.topo
        if row.p_diam_opt > 25:
            print(NAME)
        prop_file = NAME+'_opt_properties.json'
        with open(prop_file, 'r') as f:
            data = json.load(f)
        asymmetry = stk_functions.get_asymmetry(data)
        row.assym = asymmetry
        # remove structures with no. windows < expected for topology
        if row.w_no < stk_functions.expected_window(row.topo):
            continue
        # get synthetic accessibility of amines
        row.SA2 = float(SA_data[SA_data['name'] == row.bb2].SC)
        amines.append(row.bb2)
        cages.append(row['name'])
        final_db = final_db.append(row)
    return final_db


def brute_analysis(output_csv, amine_type,
                   precursor_dir, precursor_files, DB, amines,
                   macromod_, des_topo, SA_data):
    '''Function to run all analysis and screening.

    '''
    topo_markers = {'1p1': 'o', '4p4': '<', '2p3': 'X', '4p6': 'D',
                    '4p62': 'P', '6p9': '>'}
    topo_lab = {'1p1': '1+1', '4p4': '4+4', '2p3': '2+3', '4p6': '4+6',
                '4p62': 'P', '6p9': '6+9'}
    ald_colo = {'aldehyde1': 'k', 'aldehyde2': 'b',
                'aldehyde3': 'r', 'aldehyde4': 'g'}
    full_dataset = pd.read_csv(output_csv)
    working_dataset = copy.deepcopy(full_dataset)
    prec_ey_file = precursor_dir+'all_prec_ey.json'
    cage_ey_file = 'all_cage_ey.json'
    # get energies of bb1
    print('getting bb1 energies')
    energies = get_OPLS3_energy_of_list(out_file=prec_ey_file,
                                        structures=precursor_files,
                                        dir=precursor_dir,
                                        macromod_=macromod_, opt=True)
    bb1_energies = []
    for i, row in working_dataset.iterrows():
        bb1_energies.append(energies[row.bb1])
    working_dataset['bb1_ey'] = bb1_energies
    # get energies of bb2
    print('getting bb2 energies')
    energies = get_OPLS3_energy_of_list(out_file=prec_ey_file,
                                        structures=amines,
                                        dir=DB,
                                        macromod_=macromod_, opt=True)
    bb2_energies = []
    for i, row in working_dataset.iterrows():
        bb2_energies.append(energies[row.bb2])
    working_dataset['bb2_ey'] = bb2_energies
    # get energies of cages
    print('getting cage energies')
    cage_files = []
    for i, row in working_dataset.iterrows():
        NAME = row.bb1+'_'+row.bb2+'_'+row.topo
        cage_file = NAME
        cage_files.append(cage_file)
    energies = get_OPLS3_energy_of_list(out_file=cage_ey_file,
                                        structures=cage_files,
                                        dir='',
                                        macromod_=macromod_,
                                        opt=False)
    cage_energies = []
    for i, row in working_dataset.iterrows():
        NAME = row.bb1+'_'+row.bb2+'_'+row.topo
        cage_energies.append(energies[NAME])
    working_dataset['cage_ey'] = cage_energies
    working_dataset['form_ey'] = get_formation_energies(working_dataset)
    # calculate the most stable (by formation energy) of each topology
    # for a pair of BB
    # add column to final DB of relative form energy (0 if most stable topo)
    energies = []
    for i, row in working_dataset.iterrows():
        new_db = working_dataset[working_dataset.bb1 == row.bb1]
        new_db = new_db[new_db.bb2 == row.bb2]
        FEs = sorted(list(new_db.form_ey))
        energies.append(row.form_ey - min(FEs))
    working_dataset['rel_form_ey'] = energies
    # screen cages
    final_db = screening_cages(dataset=working_dataset, des_topo=des_topo,
                               SA_data=SA_data)
    print('>>> ', len(final_db), 'cages remaining after screening!')
    print('>>> doing all plotting')
    if amine_type == '2':
        plots_amines2(final_db=final_db, ald_colo=ald_colo, DB=DB,
                      des_topo=des_topo)
    if amine_type == '3':
        plots_amines3(final_db=final_db, ald_colo=ald_colo, DB=DB)
    print('>>> done!')


def main():
    """Run script.

    """
    if (not len(sys.argv) == 5):
        print("""
Usage: build_cages.py amine_type output_file wipe run_build
    amine_type: whether to use amines2f or amines3f (2 or 3).
    output_file: file to output results
    wipe: t/T if wipe output file
    run_build: t/T if you want to run the build stage""")
        sys.exit()
    else:
        amine_type = sys.argv[1]
        output_file = sys.argv[2]
        wipe = sys.argv[3]
        run_b = sys.argv[4]

    macromod_ = '/home/atarzia/software/schrodinger_install'
    precursor_dir = '/home/atarzia/projects/andrew_marsh_structures/precursor_lib/'
    # get precursor files
    precursor_files = sorted(glob.glob(precursor_dir+'*.mol'))
    precursor_names = [i.replace(precursor_dir, '') for i in precursor_files]
    # read in mol files
    precursor_struc = [stk.StructUnit3(i) for i in precursor_files]
    big_DB = '/data/atarzia/precursor_DBs/reaxys_sorted/'
    if amine_type == '2':
        amines2f = big_DB+'amines2f/'
        # synthetic accessibility DBs
        amine_SA = pd.read_csv(amines2f+'score_output_amines2f.csv')
        print(len(glob.glob(amines2f+'*mol')), 'precusrors in', amines2f)
        # diamines (alde3+amine2)
        topology_names = ['2p3', '4p6', '4p62', '6p9']  # 'dodec', '8p12',  ]
        topology_options = [stk.two_plus_three.TwoPlusThree(),
                            stk.two_plus_three.FourPlusSix(),
                            stk.two_plus_three.FourPlusSix2(),
                            stk.two_plus_three.SixPlusNine()]
                            # , stk.two_plus_three.Dodecahedron(),
                            # stk.two_plus_three.EightPlusTwelve(),]
        DB = amines2f
        amines = glob.glob(DB+'*.mol')
        prefix = 'amine2f'
        des_topo = ['2p3', '4p6', '6p9']  # '4p62',
    elif amine_type == '3':
        amines3f = big_DB+'amines3f/'
        # synthetic accessibility DBs
        amine_SA = pd.read_csv(amines3f+'score_output_amines3f.csv')
        print(len(glob.glob(amines3f+'*mol')), 'precusrors in', amines3f)
        # triamines (alde3+amine3)
        topology_names = ['1p1', '4p4']  # , '2p2']
        topology_options = [
                    stk.three_plus_three.OnePlusOne(
                        # place bb1 on vertex (0), bb2 on vertex (1)
                        bb_positions={0: [0], 1: [1]}),
                    stk.three_plus_three.FourPlusFour(
                        # place bb1 on vertex (0, 2), bb2 on vertex (1, 3)
                        bb_positions={0: [0, 3, 5, 6], 1: [1, 2, 4, 7]})
                    # stk.three_plus_three.TwoPlusTwo(
                        # place bb1 on vertex (0, 2), bb2 on vertex (1, 3)
                        # bb_positions={0: [0, 2], 1: [1, 3]})]
                    ]
        DB = amines3f
        amines = glob.glob(DB+'*.mol')
        prefix = 'amine3f'
        des_topo = ['1p1', '4p4']

    # prepare output file
    output_csv = output_file
    if wipe.lower() == 't':
        with open(output_csv, 'w') as f:
            f.write('name,bb1,SA1,bb2,SA2,topo,max_diam,p_diam,p_vol,p_diam_opt')
            f.write(',p_vol_opt,w_no,w_max,w_min,w_avg,w_diff,collapse,assym\n')
    if run_b.lower() == 't':
        brute_cage_build(precursor_struc, precursor_names, precursor_files,
                         topology_names, topology_options, amines, DB,
                         amine_type, output_csv, macromod_)
    check_output_numbers(precursor_dir, DB, topology_options, prefix)
    brute_analysis(output_csv, amine_type,
                   precursor_dir, precursor_files, DB, amines,
                   macromod_, des_topo, amine_SA)


if __name__ == "__main__":
    main()
