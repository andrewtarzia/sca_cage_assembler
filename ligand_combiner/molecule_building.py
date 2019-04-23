#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build all molecules from DB of core, ligands and linkers.

Author: Andrew Tarzia

Date Created: 04 Apr 2019

"""

import sys
import pickle
import numpy as np
from itertools import product
import pandas as pd
from rdkit.Chem import AllChem as Chem
from os.path import join, isfile
import stk
from stk import Population
from Combiner import get_molecule, get_geometrical_properties, Combination
sys.path.insert(0, '/home/atarzia/thesource/')
from stk_functions import build_population
from rdkit_functions import mol_list2grid
from analyze_molecules import plot_all_pair_info


def main():
    if (not len(sys.argv) == 8):
        print("""
    Usage: molecule_builing.py rebuild N bond_mean bond_std pair_data energy_tol angle_tol
        rebuild (str) - 't' if you want to rebuild all molecules, 'f' to load populations
        N (int) - number of conformers to use
        bond_mean (float) - mean value of bond distance to use in candidate selection
        bond_std (float) - std deviation value of bond distance to use in candidate selection
        pair_data (str) - pickle file to output pair data to
        energy_tol (float) - max kJ/mol over min energy conformer to allow
        angle_tol (float) - tolerance to use on angle matching
        """)
        sys.exit()
    else:
        rebuild = sys.argv[1]
        N = int(sys.argv[2])
        bond_mean = float(sys.argv[3])
        bond_std = float(sys.argv[4])
        pair_data = str(sys.argv[5])
        energy_tol = float(sys.argv[6])
        angle_tol = float(sys.argv[7])

    proj_dir = '/home/atarzia/projects/ligand_combiner/'
    core_dir = proj_dir + 'cores/'
    liga_dir = proj_dir + 'ligands/'
    link_dir = proj_dir + 'linkers/'
    mole_dir = proj_dir + 'molecules/'

    if rebuild == 't':
        print('building molecules')
        # build molecule populations
        core_pop = build_population(directory=core_dir, structunit='StructUnit2',
                                    fgs=['bromine'])
        liga_pop = build_population(directory=liga_dir, structunit='StructUnit',
                                    fgs=['bromine'])
        link_pop = build_population(directory=link_dir, structunit='StructUnit2',
                                    fgs=['bromine'])

        # for the linker molecules, we attach an invertable flag attribute
        link_mol_prop = pd.read_csv(join(link_dir, 'linkers.csv'))
        for i in link_pop:
            inversion_flag = str(link_mol_prop[link_mol_prop.name == i.name]['inversion_flag'].iloc[0])
            if inversion_flag == 't':
                i.invertable = True
            else:
                i.invertable = False

        # this is the resultant molecule population
        molecule_pop = Population()
        for item in product(enumerate(core_pop), enumerate(liga_pop), enumerate(link_pop)):
            core_item, liga_item, link_item = item
            i, core = core_item
            j, liga = liga_item
            k, link = link_item
            if core.name not in ['core_4', 'core_5', 'core_6']:
                # if core.name not in ['core_4', 'core_6']:
                continue
            if liga.name not in ['lig_1', 'lig_2', 'lig_3']:
                # if liga.name not in ['lig_1', 'lig_2']:
                continue
            if link.name not in ['link_1']:
                continue
            # build ABCBA molecule
            pop_ids = (i, j, k)  # core, ligand, linker
            print(core.name, liga.name, link.name, pop_ids)
            ABCBA_confs, ABCBA_molecule = get_molecule(
                type='ABCBA', inverted=False,
                popns=(core_pop, liga_pop, link_pop),
                pop_ids=pop_ids, N=N,
                mole_dir=mole_dir)
            ABCBA_molecule.name = 'ABCBA_' + str(i) + str(j) + str(k)
            # get properties - save to molecule as attribute
            ABCBA_molecule = get_geometrical_properties(mol=ABCBA_molecule,
                                                        cids=ABCBA_confs,
                                                        type='ABCBA')
            molecule_pop.members.append(ABCBA_molecule)
            # also build the inverted molecule if possible.
            if link_pop[pop_ids[2]].invertable is True:
                ABCBA_inv_confs, ABCBA_inv_molecule = get_molecule(
                    type='ABCBA', inverted=True,
                    popns=(core_pop, liga_pop, link_pop),
                    pop_ids=pop_ids, N=N,
                    mole_dir=mole_dir)
                ABCBA_inv_molecule.name = 'ABCBA_' + str(i) + str(j) + str(k) + 'i'
                # get properties - save to molecule as attribute
                ABCBA_inv_molecule = get_geometrical_properties(mol=ABCBA_inv_molecule,
                                                                cids=ABCBA_inv_confs,
                                                                type='ABCBA')
                molecule_pop.members.append(ABCBA_inv_molecule)
            # avoid doing multple times
            if isfile(core.name + '_' + liga.name + '_opt.mol') is False:
                # build ABA molecule
                ABA_confs, ABA_molecule = get_molecule(
                    type='ABA', inverted=False,
                    popns=(core_pop, liga_pop, link_pop),
                    pop_ids=pop_ids, N=N,
                    mole_dir=mole_dir)
                ABA_molecule.name = 'ABA_' + str(i) + str(j)
                # get properties - save to molecule as attribute
                ABA_molecule = get_geometrical_properties(mol=ABA_molecule,
                                                          cids=ABA_confs,
                                                          type='ABA')
                molecule_pop.members.append(ABA_molecule)

        # save populations
        core_pop.dump(join(core_dir, 'core.pop'), include_attrs=['geom_prop'])
        liga_pop.dump(join(liga_dir, 'ligands.pop'), include_attrs=['geom_prop'])
        link_pop.dump(join(link_dir, 'linkers.pop'), include_attrs=['geom_prop'])
        molecule_pop.dump(join(mole_dir, 'molecules.pop'), include_attrs=['geom_prop'])

        # draw 2D representation for all built molecules
        mol_list = []
        for poly in molecule_pop:
            MOL = Chem.MolFromSmiles(Chem.MolToSmiles(poly.mol))
            mol_list.append(MOL)
        mol_list2grid(molecules=mol_list, filename='built_molecules',
                      mol_per_row=3, maxrows=3, subImgSize=(200, 200))
        print('done')
        print('----------------------------------')
    elif rebuild == 'f':
        print('loading in populations')
        # load in populations
        core_pop = stk.Population.load(path=join(core_dir, 'core.pop'),
                                       member_init=stk.Molecule.from_dict)
        liga_pop = stk.Population.load(path=join(liga_dir, 'ligands.pop'),
                                       member_init=stk.Molecule.from_dict)
        link_pop = stk.Population.load(path=join(link_dir, 'linkers.pop'),
                                       member_init=stk.Molecule.from_dict)
        molecule_pop = stk.Population.load(path=join(mole_dir, 'molecules.pop'),
                                           member_init=stk.Molecule.from_dict)
        print('done')
        print('----------------------------------')

    print('obtaining properties for all pairs')
    # define bond length vector to use based on N-Pd bond distances extracted
    # from survey
    # N-Pd-N length
    vector_length = 2 * bond_mean
    vector_std = 2 * bond_std

    # obtain all pair properties in molecule DB
    # poly1 should be the 'large' molecule, while poly2 should be the 'small'
    # molecule of the pair
    all_pairs = []
    for i, poly1 in enumerate(molecule_pop):
        print('i', i)
        for j, poly2 in enumerate(molecule_pop):
            #################################
            # for setting specific polymers
            # if poly2.name != 'ABA_51':
            #     continue
            # if poly1.name != 'ABCBA_300':
            #     continue
            #################################
            # make sure poly1 != poly2
            if i == j:
                continue
            for conf1 in poly1.geom_prop:
                PROP1 = poly1.geom_prop[conf1]
                # skip conformer if dihedral meant the N's were not on the
                # right side or Ns are pointing the wrong way
                if PROP1['skip'] is True:
                    continue
                for conf2 in poly2.geom_prop:
                    PROP2 = poly2.geom_prop[conf2]
                    # skip conformer if dihedral meant the N's were not on the
                    # right side or Ns are pointing the wrong way
                    if PROP2['skip'] is True:
                        continue
                    comb = Combination(poly1, poly2, conf1, conf2)
                    comb.popn_ids = (i, j)
                    # if molecule1 or molecule2 energy > threshold from conf min
                    # skip pair
                    if comb.energy1 > energy_tol or comb.energy2 > energy_tol:
                        print('skipping conformer pair due to energy')
                        print(i, conf1, j, conf2)
                        print(comb.energy1, comb.energy2)
                        continue
                    # obtain all properties
                    # check N-N distance of poly1-conf > poly2-conf
                    comb.NN_dist1 = np.linalg.norm(np.asarray(PROP1['NN_v']))
                    comb.NN_dist2 = np.linalg.norm(np.asarray(PROP2['NN_v']))
                    # only save combinations with NN_dist1 > NN_dist2
                    if comb.test_N_N_lengths is False:
                        continue
                    # determine angles made by NN_v and NN-BC_v
                    # check that the pairs sum to 180
                    comb.p1_angle1 = PROP1['NN_BCN_1']
                    comb.p1_angle2 = PROP1['NN_BCN_2']
                    comb.p2_angle1 = PROP2['NN_BCN_1']
                    comb.p2_angle2 = PROP2['NN_BCN_2']
                    # now check that the length of the long
                    # vector and the short vector are commensurate
                    # with an ideal trapezoid with the given angles
                    # i.e. the extender vector determined by the
                    # difference of the two NN_v (LHS) matches what is
                    # expected by trig (RHS)
                    comb.extender_V_LHS = (comb.NN_dist1 - comb.NN_dist2) / 2
                    comb.test_angle = np.radians(180 - comb.p2_angle1)
                    comb.extender_V_RHS = vector_length * np.cos(comb.test_angle)
                    comb.tol = vector_std * np.cos(comb.test_angle)
                    # get final geometrical properties
                    comb.get_angle_deviations()
                    comb.get_N_Pd_lengths_deviation()
                    all_pairs.append(comb)
    print('done')
    print(len(all_pairs))
    print('----------------------------------')
    # do analysis
    print('doing all analysis')
    plot_all_pair_info(pair_data=all_pairs,
                       angle_tol=angle_tol, energy_tol=energy_tol,
                       outfile=pair_data.replace('.pkl', '.pdf'))
    print('done')
    print('----------------------------------')


if __name__ == "__main__":
    main()
