#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build all molecules from DB of core, ligands and linkers.

Author: Andrew Tarzia

Date Created: 04 Apr 2019

"""

import sys
import itertools
import pandas as pd
from rdkit.Chem import AllChem as Chem
import os
import glob
import stk
sys.path.insert(0, '/home/atarzia/thesource/')
import Combiner
import stk_functions
import rdkit_functions


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
    if (not len(sys.argv) == 4):
        print("""
    Usage: molecule_builing.py subset rebuild N bond_mean bond_std energy_tol angle_tol
        subset (str) - 'all' to build all possible molecules,
            'clever' to build molecules from bloch2017
        rebuild (str) - 't' if you want to rebuild all molecules, 'f' to load populations
        N (int) - number of conformers to use
        """)
        sys.exit()
    else:
        subset = sys.argv[1]
        rebuild = sys.argv[2]
        N = int(sys.argv[3])

    proj_dir = '/home/atarzia/projects/ligand_combiner/'
    core_dir = proj_dir + 'cores/'
    liga_dir = proj_dir + 'ligands/'
    link_dir = proj_dir + 'linkers/'

    if rebuild == 't':
        print('building molecules')
        # clear already built molecules
        for file in glob.glob('AB*.pdb'):
            os.remove(os.path.join('./', file))
        for file in glob.glob('AB*_POIs.xyz'):
            os.remove(os.path.join('./', file))
        for file in glob.glob('core*.mol'):
            os.remove(os.path.join('./', file))
        for file in glob.glob('core*.json'):
            os.remove(os.path.join('./', file))
        for file in glob.glob('built_molecules*'):
            os.remove(os.path.join('./', file))
        # build molecule populations
        core_pop = stk_functions.build_population(directory=core_dir,
                                                  structunit='StructUnit2',
                                                  fgs=['bromine'])
        liga_pop = stk_functions.build_population(directory=liga_dir,
                                                  structunit='StructUnit',
                                                  fgs=['bromine'])
        link_pop = stk_functions.build_population(directory=link_dir,
                                                  structunit='StructUnit2',
                                                  fgs=['bromine'])

        # for the linker molecules, we attach an invertable flag attribute
        link_mol_prop = pd.read_csv(os.path.join(link_dir, 'linkers.csv'))
        for i in link_pop:
            inversion_flag = str(link_mol_prop[link_mol_prop.name == i.name]['inversion_flag'].iloc[0])
            if inversion_flag == 't':
                i.invertable = True
            else:
                i.invertable = False

        # this is the resultant molecule population
        molecule_pop = stk.Population()
        iterator = itertools.product(enumerate(core_pop),
                                     enumerate(liga_pop),
                                     enumerate(link_pop))
        for item in iterator:
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
                ABCBA_confs, ABCBA_molecule = Combiner.get_molecule(
                    type='ABCBA', inverted=False,
                    popns=(core_pop, liga_pop, link_pop),
                    pop_ids=pop_ids, N=N)
                ABCBA_molecule.name = 'ABCBA_' + str(i) + str(j) + str(k)
                # get properties - save to molecule as attribute
                ABCBA_molecule = Combiner.get_geometrical_properties(mol=ABCBA_molecule,
                                                                     cids=ABCBA_confs,
                                                                     type='ABCBA')
                molecule_pop.members.append(ABCBA_molecule)
            # also build the inverted molecule if possible.
            # if building only a subset, check that this molecule is in the subset
            NAME = core.name + '_' + liga.name + '_' + link.name + 'i'
            if subset != 'all' and NAME in cases(subset):
                if link_pop[pop_ids[2]].invertable is True:
                    ABCBA_inv_confs, ABCBA_inv_molecule = Combiner.get_molecule(
                        type='ABCBA', inverted=True,
                        popns=(core_pop, liga_pop, link_pop),
                        pop_ids=pop_ids, N=N)
                    ABCBA_inv_molecule.name = 'ABCBA_' + str(i) + str(j) + str(k) + 'i'
                    # get properties - save to molecule as attribute
                    ABCBA_inv_molecule = Combiner.get_geometrical_properties(mol=ABCBA_inv_molecule,
                                                                             cids=ABCBA_inv_confs,
                                                                             type='ABCBA')
                    molecule_pop.members.append(ABCBA_inv_molecule)
            # build ABA molecule
            # avoid doing multple times
            # if building only a subset, check that this molecule is in the subset
            NAME = core.name + '_' + liga.name
            if subset != 'all' and NAME in cases(subset):
                if os.path.isfile(core.name+'_'+liga.name+'_ABA_opt.mol') is False:
                    print(os.path.isfile(core.name+'_'+liga.name+'_ABA_opt.mol'), core.name+'_'+liga.name+'_ABA_opt.mol')
                    ABA_confs, ABA_molecule = Combiner.get_molecule(
                        type='ABA', inverted=False,
                        popns=(core_pop, liga_pop, link_pop),
                        pop_ids=pop_ids, N=N)
                    ABA_molecule.name = 'ABA_' + str(i) + str(j)
                    # get properties - save to molecule as attribute
                    ABA_molecule = Combiner.get_geometrical_properties(mol=ABA_molecule,
                                                                       cids=ABA_confs,
                                                                       type='ABA')
                    molecule_pop.members.append(ABA_molecule)

        # save populations
        core_pop.dump(os.path.join(core_dir, 'core.pop'), include_attrs=['geom_prop'])
        liga_pop.dump(os.path.join(liga_dir, 'ligands.pop'), include_attrs=['geom_prop'])
        link_pop.dump(os.path.join(link_dir, 'linkers.pop'), include_attrs=['geom_prop'])
        molecule_pop.dump(os.path.join('./', 'molecules.pop'), include_attrs=['geom_prop'])

        # draw 2D representation for all built molecules
        mol_list = []
        for poly in molecule_pop:
            MOL = Chem.MolFromSmiles(Chem.MolToSmiles(poly.mol))
            mol_list.append(MOL)
        rdkit_functions.mol_list2grid(molecules=mol_list, filename='built_molecules',
                                      mol_per_row=3, maxrows=3, subImgSize=(200, 200))
        print('done')
        print('----------------------------------')


if __name__ == "__main__":
    main()
