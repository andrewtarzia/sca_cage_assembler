#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build all molecules from DB of core, ligands and linkers.

Author: Andrew Tarzia

Date Created: 04 Apr 2019

"""

import sys
from itertools import product
from pandas import read_csv
from rdkit.Chem import AllChem as Chem
from os.path import join, isfile
from os import remove
from glob import glob
from stk import Population
from Combiner import get_molecule, get_geometrical_properties
sys.path.insert(0, '/home/atarzia/thesource/')
from stk_functions import build_population
from rdkit_functions import mol_list2grid


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


if __name__ == "__main__":
    main()
