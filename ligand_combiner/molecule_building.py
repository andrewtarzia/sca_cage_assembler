#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build all molecules from DB of core, ligands and linkers.

Author: Andrew Tarzia

Date Created: 04 Apr 2019

"""

import sys
from ase.atoms import Atoms, Atom
import numpy as np
import pandas as pd
from rdkit.Chem import AllChem as Chem
from os.path import join
import matplotlib.pyplot as plt
import stk
from stk import rdkit_ETKDG, Population
sys.path.insert(0, '/home/atarzia/thesource/')
from stk_functions import build_population, build_ABCBA, build_ABA
from calculations import get_dihedral, angle_between
from rdkit_functions import mol_list2grid


def main():
    if (not len(sys.argv) == 6):
        print("""
    Usage: molecule_builing.py rebuild N angle_tol bond_mean bond_std
        rebuild (str) - 't' if you want to rebuild all molecules, 'f' to load populations
        N (int) - number of conformers to use
        angle_tol (float) - tolerance to use on angle matching
        bond_mean (float) - mean value of bond distance to use in candidate selection
        bond_std (float) - std deviation value of bond distance to use in candidate selection
        """)
        sys.exit()
    else:
        rebuild = sys.argv[1]
        N = int(sys.argv[2])
        angle_tol = float(sys.argv[3])
        bond_mean = float(sys.argv[4])
        bond_std = float(sys.argv[5])

    proj_dir = '/home/atarzia/projects/ligand_combiner/'
    core_dir = proj_dir + 'cores/'
    liga_dir = proj_dir + 'ligands/'
    link_dir = proj_dir + 'linkers/'
    mole_dir = proj_dir + 'molecules/'

    if rebuild == 't':
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
        for i, core in enumerate(core_pop):
            # if core.name not in ['core_4', 'core_5', 'core_6']:
            if core.name not in ['core_4', 'core_6']:
                continue
            for j, liga in enumerate(liga_pop):
                # if liga.name not in ['lig_1', 'lig_2', 'lig_3']:
                if liga.name not in ['lig_1', 'lig_2']:
                    continue
                for k, link in enumerate(link_pop):
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
                    # break
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
                # break
            # break

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

    elif rebuild == 'f':
        # load in populations
        core_pop = stk.Population.load(path=join(core_dir, 'core.pop'),
                                       member_init=stk.Molecule.from_dict)
        liga_pop = stk.Population.load(path=join(liga_dir, 'ligands.pop'),
                                       member_init=stk.Molecule.from_dict)
        link_pop = stk.Population.load(path=join(link_dir, 'linkers.pop'),
                                       member_init=stk.Molecule.from_dict)
        molecule_pop = stk.Population.load(path=join(mole_dir, 'molecules.pop'),
                                           member_init=stk.Molecule.from_dict)

    # define bond length vector to use based on N-Pd bond distances extracted
    # from survey
    # N-Pd-N length
    vector_length = 2 * bond_mean
    vector_std = 2 * bond_std

    # find matching pairs in molecule DB
    # poly1 should be the 'large' molecule, while poly2 should be the 'small'
    # molecule of the pair
    # this is a crude/lazy for loop
    passed_pairs = []
    all_pairs = []
    for i, poly1 in enumerate(molecule_pop):
        for conf1 in poly1.geom_prop:
            if poly1.name != 'ABCBA_300':
                continue
            PROP1 = poly1.geom_prop[conf1]
            # skip conformer if dihedral meant the N's were not on the right side
            if PROP1['skip'] is True:
                continue
            print('doing conformer {} of ABCBA'.format(conf1))
            for j, poly2 in enumerate(molecule_pop):
                if poly2.name != 'ABA_51':
                    continue
                # make sure poly1 != poly2
                if i != j:
                    NN_dist2s = []
                    angles1_sums = []
                    angles2_sums = []
                    LHSs = []
                    RHSs = []
                    for conf2 in poly2.geom_prop:
                        PROP2 = poly2.geom_prop[conf2]
                        # skip conformer if desired
                        if PROP2['skip'] is True:
                            continue
                        print(poly1.name, conf1, poly2.name, conf2)
                        # obtain all properties
                        # check N-N distance of poly1-conf > poly2-conf
                        NN_dist1 = np.linalg.norm(PROP1['NN_v'])
                        NN_dist2 = np.linalg.norm(PROP2['NN_v'])
                        # determine angles made by NN_v and NN-BC_v
                        # check that the pairs sum to 180
                        p1_angle1 = PROP1['NN_BCN_1']
                        p1_angle2 = PROP1['NN_BCN_2']
                        p2_angle1 = PROP2['NN_BCN_1']
                        p2_angle2 = PROP2['NN_BCN_2']
                        # now check that the length of the long
                        # vector and the short vector are commensurate
                        # with an ideal trapezoid with the given angles
                        # i.e. the extender vector determined by the
                        # difference of the two NN_v (LHS) matches what is
                        # expected by trig (RHS)
                        extender_V_LHS = (NN_dist1 - NN_dist2) / 2
                        test_angle = np.radians(180 - p2_angle1)
                        extender_V_RHS = vector_length * np.cos(test_angle)
                        tol = 2 * vector_std * np.cos(test_angle)
                        # temporary for plotting
                        NN_dist2s.append(NN_dist2)
                        angles1_sums.append(p1_angle1 + p2_angle1)
                        angles2_sums.append(p1_angle2 + p2_angle2)
                        LHSs.append(extender_V_LHS)
                        RHSs.append(extender_V_RHS)
                        # run checks
                        if NN_dist1 > NN_dist2:
                            if np.isclose(p1_angle1 + p2_angle1, 180, rtol=0, atol=angle_tol):
                                if np.isclose(p1_angle2 + p2_angle2, 180, rtol=0, atol=angle_tol):
                                    if np.isclose(extender_V_LHS, extender_V_RHS,
                                                  rtol=0, atol=tol):
                                        print(poly1.name, conf1, poly2.name, conf2)
                                        print('NN_dists:', NN_dist1, NN_dist2)
                                        print('angles:', p1_angle1,
                                              p1_angle2, p2_angle1, p2_angle2)
                                        print('sum angles:',
                                              p1_angle1 + p2_angle1,
                                              p1_angle2 + p2_angle2)
                                        print('LHS:', extender_V_LHS)
                                        print('cos(pi-d):', test_angle)
                                        print('RHS:', extender_V_RHS)
                                        print('bond_tol:', tol)
                                        print('passed')
                                        input()
                                        passed_pairs.append((i, j, conf1, conf2,
                                                             p1_angle1, p1_angle2,
                                                             p2_angle1, p2_angle2,
                                                             NN_dist1, NN_dist2))
                            all_pairs.append((i, j, conf1, conf2,
                                              p1_angle1, p1_angle2,
                                              p2_angle1, p2_angle2,
                                              NN_dist1, NN_dist2))
                    # make some temp plots
                    fig, ax = plt.subplots(figsize=(8, 5))
                    X_range = (5, 20)
                    width = 0.5
                    density = False
                    Y = NN_dist2s
                    xtitle = 'NN_dist2'
                    ax.axvline(x=NN_dist1, c='k')
                    X_bins = np.arange(X_range[0], X_range[1], width)
                    hist, bin_edges = np.histogram(a=Y, bins=X_bins,
                                                   density=density)
                    ax.bar(bin_edges[:-1],
                           hist,
                           align='edge',
                           alpha=1.0, width=width,
                           color='firebrick',
                           edgecolor='k')
                    ax.tick_params(axis='both', which='major', labelsize=16)
                    ax.set_xlabel(xtitle, fontsize=16)
                    if density is False:
                        ax.set_ylabel('count', fontsize=16)
                    elif density is True:
                        ax.set_ylabel('frequency', fontsize=16)
                    ax.set_xlim(X_range)
                    ax.set_title(poly1.name + '_' + str(conf1))
                    fig.tight_layout()
                    fig.savefig('NN_dists_' + poly1.name + '_' + str(conf1) + '_'+ poly2.name + '.pdf', dpi=720,
                                bbox_inches='tight')
                    plt.close()

                    fig, ax = plt.subplots(figsize=(8, 5))
                    X_range = (90, 270)
                    width = 5
                    density = False
                    xtitle = 'sum of angles [degrees]'
                    ax.axvline(x=180 - angle_tol, c='k')
                    ax.axvline(x=180 + angle_tol, c='k')

                    Y = angles1_sums
                    X_bins = np.arange(X_range[0], X_range[1], width)
                    hist, bin_edges = np.histogram(a=Y, bins=X_bins,
                                                   density=density)
                    ax.bar(bin_edges[:-1],
                           hist,
                           align='edge',
                           alpha=0.6, width=width,
                           color='r',
                           edgecolor='k', label='angle1 sum')
                    Y = angles2_sums
                    X_bins = np.arange(X_range[0], X_range[1], width)
                    hist, bin_edges = np.histogram(a=Y, bins=X_bins,
                                                   density=density)
                    ax.bar(bin_edges[:-1],
                           hist,
                           align='edge',
                           alpha=0.6, width=width,
                           color='b',
                           edgecolor='k', label='angle2 sum')

                    ax.tick_params(axis='both', which='major', labelsize=16)
                    ax.set_xlabel(xtitle, fontsize=16)
                    if density is False:
                        ax.set_ylabel('count', fontsize=16)
                    elif density is True:
                        ax.set_ylabel('frequency', fontsize=16)
                    ax.set_xlim(X_range)
                    ax.set_title(poly1.name + '_' + str(conf1))
                    ax.legend()
                    fig.tight_layout()
                    fig.savefig('angles_dists_' + poly1.name + '_' + str(conf1) + '_'+ poly2.name + '.pdf', dpi=720,
                                bbox_inches='tight')
                    plt.close()

                    fig, ax = plt.subplots(figsize=(5, 5))
                    X_range = (0, 5)
                    X = LHSs
                    Y = RHSs
                    ax.scatter(X, Y, alpha=0.6, color='r', edgecolor='k')
                    ax.plot(np.linspace(0, X_range[1], 10),
                            np.linspace(0, X_range[1], 10), alpha=0.2, c='k')
                    ax.tick_params(axis='both', which='major', labelsize=16)
                    ax.set_xlabel('LHS', fontsize=16)
                    ax.set_ylabel('RHS', fontsize=16)
                    ax.set_xlim(X_range)
                    ax.set_ylim(X_range)
                    ax.set_title(poly1.name + '_' + str(conf1))
                    fig.tight_layout()
                    fig.savefig('LHSRHS_parity_' + poly1.name + '_' + str(conf1) + '_'+ poly2.name + '.pdf', dpi=720,
                                bbox_inches='tight')
                    plt.close()


if __name__ == "__main__":
    main()
