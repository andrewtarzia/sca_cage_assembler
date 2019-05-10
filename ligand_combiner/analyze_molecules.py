#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build all molecules from DB of core, ligands and linkers.

Author: Andrew Tarzia

Date Created: 05 May 2019

"""

import sys
from numpy.linalg import norm
from numpy import asarray, cos, radians
from itertools import product
from pandas import read_csv
from rdkit.Chem import AllChem as Chem
from os.path import join, isfile
from os import remove
from glob import glob
from stk import Population, Molecule
from Combiner import get_molecule, get_geometrical_properties, Combination
sys.path.insert(0, '/home/atarzia/thesource/')
from stk_f import build_population
from rdkit_f import mol_list2grid
from analysis import plot_all_pair_info, output_analysis


def cases(subset):
    '''Define the smallest set of molecules to build to replicate particular subsets

    Implemented:
        bloch2017:
            1 (1 in Fig 4) — (LA) core = 5, linker = 1, ligand = 3 + (LP) core = 6, ligand = 2
            2 (2 in Fig 4) — (LC) core = 4, linker = 1, ligand = 1 + (LP) core = 6, ligand = 2
            3 (3 in Fig 4) — (LA) core = 5, linker = 1, ligand = 3  + (LC) core = 4, linker = 1, ligand = 1 — named interlocked_*

    '''
    cases = {'clever':
             ['core_5_lig_3_link_1', 'core_6_lig_2', 'core_4_lig_1_link_1']
             }
    try:
        list_of_mol = cases[subset]
    except KeyError:
        print('{} is not a defined subset'.format(subset))
        sys.exit('exitting')
    return list_of_mol


def main():
    if (not len(sys.argv) == 8):
        print("""
    Usage: molecule_builing.py subset rebuild N bond_mean bond_std energy_tol angle_tol
        subset (str) - 'all' to build all possible molecules,
            'clever' to build molecules from bloch2017
        rebuild (str) - 't' if you want to rebuild all molecules, 'f' to load populations
        N (int) - number of conformers to use
        bond_mean (float) - mean value of bond distance to use in candidate selection
        bond_std (float) - std deviation value of bond distance to use in candidate selection
        energy_tol (float) - max kJ/mol over min energy conformer to allow
        angle_tol (float) - tolerance to use on angle matching
        """)
        sys.exit()
    else:
        subset = sys.argv[1]
        rebuild = sys.argv[2]
        N = int(sys.argv[3])
        bond_mean = float(sys.argv[4])
        bond_std = float(sys.argv[5])
        energy_tol = float(sys.argv[6])
        angle_tol = float(sys.argv[7])

    proj_dir = '/home/atarzia/projects/ligand_combiner/'
    core_dir = proj_dir + 'cores/'
    liga_dir = proj_dir + 'ligands/'
    link_dir = proj_dir + 'linkers/'

    if rebuild == 't':
        print('building molecules')
        # clear already built molecules
        for file in glob('AB*.pdb'):
            remove(join('./', file))
        for file in glob('AB*_POIs.xyz'):
            remove(join('./', file))
        for file in glob('core*.mol'):
            remove(join('./', file))
        for file in glob('core*.json'):
            remove(join('./', file))
        for file in glob('built_molecules*'):
            remove(join('./', file))
        # build molecule populations
        core_pop = build_population(directory=core_dir, structunit='StructUnit2',
                                    fgs=['bromine'])
        liga_pop = build_population(directory=liga_dir, structunit='StructUnit',
                                    fgs=['bromine'])
        link_pop = build_population(directory=link_dir, structunit='StructUnit2',
                                    fgs=['bromine'])

        # for the linker molecules, we attach an invertable flag attribute
        link_mol_prop = read_csv(join(link_dir, 'linkers.csv'))
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
            pop_ids = (i, j, k)  # core, ligand, linker
            # build ABCBA molecule
            # if building only a subset, check that this molecule is in the subset
            NAME = core.name + '_' + liga.name + '_' + link.name
            if subset != 'all' and NAME in cases(subset):
                print(core.name, liga.name, link.name, pop_ids)
                ABCBA_confs, ABCBA_molecule = get_molecule(
                    type='ABCBA', inverted=False,
                    popns=(core_pop, liga_pop, link_pop),
                    pop_ids=pop_ids, N=N)
                ABCBA_molecule.name = 'ABCBA_' + str(i) + str(j) + str(k)
                # get properties - save to molecule as attribute
                ABCBA_molecule = get_geometrical_properties(mol=ABCBA_molecule,
                                                            cids=ABCBA_confs,
                                                            type='ABCBA')
                molecule_pop.members.append(ABCBA_molecule)
            # also build the inverted molecule if possible.
            # if building only a subset, check that this molecule is in the subset
            NAME = core.name + '_' + liga.name + '_' + link.name + 'i'
            if subset != 'all' and NAME in cases(subset):
                if link_pop[pop_ids[2]].invertable is True:
                    ABCBA_inv_confs, ABCBA_inv_molecule = get_molecule(
                        type='ABCBA', inverted=True,
                        popns=(core_pop, liga_pop, link_pop),
                        pop_ids=pop_ids, N=N)
                    ABCBA_inv_molecule.name = 'ABCBA_' + str(i) + str(j) + str(k) + 'i'
                    # get properties - save to molecule as attribute
                    ABCBA_inv_molecule = get_geometrical_properties(mol=ABCBA_inv_molecule,
                                                                    cids=ABCBA_inv_confs,
                                                                    type='ABCBA')
                    molecule_pop.members.append(ABCBA_inv_molecule)
            # build ABA molecule
            # avoid doing multple times
            # if building only a subset, check that this molecule is in the subset
            NAME = core.name + '_' + liga.name
            if subset != 'all' and NAME in cases(subset):
                if isfile(core.name+'_'+liga.name+'_ABA_opt.mol') is False:
                    print(isfile(core.name+'_'+liga.name+'_ABA_opt.mol'), core.name+'_'+liga.name+'_ABA_opt.mol')
                    ABA_confs, ABA_molecule = get_molecule(
                        type='ABA', inverted=False,
                        popns=(core_pop, liga_pop, link_pop),
                        pop_ids=pop_ids, N=N)
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
        molecule_pop.dump(join('./', 'molecules.pop'), include_attrs=['geom_prop'])

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
        core_pop = Population.load(path=join(core_dir, 'core.pop'),
                                   member_init=Molecule.from_dict)
        liga_pop = Population.load(path=join(liga_dir, 'ligands.pop'),
                                   member_init=Molecule.from_dict)
        link_pop = Population.load(path=join(link_dir, 'linkers.pop'),
                                   member_init=Molecule.from_dict)
        molecule_pop = Population.load(path=join('./', 'molecules.pop'),
                                       member_init=Molecule.from_dict)
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
        print('molecule', poly1.name)
        for j, poly2 in enumerate(molecule_pop):
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
                        # print('skipping conformer pair due to energy')
                        # print(poly1.name, conf1, poly2.name, conf2)
                        # print(comb.energy1, comb.energy2)
                        continue
                    # obtain all properties
                    # check N-N distance of poly1-conf > poly2-conf
                    comb.NN_dist1 = norm(asarray(PROP1['NN_v']))
                    comb.NN_dist2 = norm(asarray(PROP2['NN_v']))
                    # only save combinations with NN_dist1 > NN_dist2
                    if comb.test_N_N_lengths() is False:
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
                    comb.test_angle = radians(180 - comb.p2_angle1)
                    comb.extender_V_RHS = vector_length * cos(comb.test_angle)
                    comb.tol = vector_std * cos(comb.test_angle)
                    # get final geometrical properties
                    comb.get_angle_deviations()
                    comb.get_N_Pd_lengths_deviation()
                    all_pairs.append(comb)
    print('done')
    print(len(all_pairs))
    print('----------------------------------')
    # do analysis
    print('doing all analysis')
    # plot_all_pair_info(pair_data=all_pairs,
    #                    angle_tol=angle_tol, energy_tol=energy_tol)
    output_analysis(molecule_pop=molecule_pop, pair_data=all_pairs,
                    angle_tol=angle_tol, energy_tol=energy_tol)
    print('done')
    print('----------------------------------')



###################
    for p, poly1 in enumerate(molecule_pop):
        combinations = [i for i in pair_data
                        if p == i.popn_ids[0]]
        print(p, poly1.name)
        print(len(combinations))
        if len(combinations) == 0:
            continue
        prefix = poly1.name + '_analysis_'

        print(combinations[0].popn_ids)
        X_data = [max([i.angle1_deviation, i.angle2_deviation])
                  for i in combinations if i.test_N_N_lengths]
        Y_data = [i.NPdN_difference
                  for i in combinations if i.test_N_N_lengths]
        Z_data = [max([i.energy1, i.energy2])/energy_tol
                  for i in combinations if i.test_N_N_lengths]
        # define colour map based on energy tol
        cmap = {'mid_point': energy_tol/2/energy_tol,
                'cmap': RdBu,
                'ticks': [0, energy_tol/2/energy_tol, energy_tol/energy_tol],
                'labels': ['0', str(energy_tol/2), str(energy_tol)],
                'cmap_label': 'energy [kJ/mol]'}

        scatter_plot(X=X_data, Y=Y_data,
                     outfile=prefix+'main.pdf',
                     xtitle='maximum angle deviation [degrees]',
                     ytitle='deviation from planarity',
                     title=poly1.name,
                     xlim=(0, 180), ylim=(0, round(max(Y_data))+1),
                     c='firebrick', edgecolors='k',
                     marker='o', alpha=0.5, s=80, Z=Z_data,
                     cmap=cmap)


if __name__ == "__main__":
    main()
