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
from stk import rdkit_ETKDG, Population
sys.path.insert(0, '/home/atarzia/thesource/')
from stk_functions import build_population, build_ABCBA, build_ABA
from calculations import get_dihedral, angle_between
from rdkit_functions import mol_list2grid


def visualize_atoms(mol, conf_dict, cid, type):
    '''Output files with POIs and conformer atom positions for visulization.

    '''
    if type == 'ABCBA':
        mol.write(path='testing_' + str(cid) + '_mol.pdb', conformer=cid)
        POIs = Atoms()
        POIs.append(Atom(symbol='H', position=conf_dict['COM']))
        POIs.append(Atom(symbol='C', position=conf_dict['liga1']['pos']))
        # POIs.append(Atom(symbol='C', position=conf_dict['link1']['pos']))
        POIs.append(Atom(symbol='C', position=conf_dict['core1']['pos']))
        # POIs.append(Atom(symbol='C', position=conf_dict['link2']['pos']))
        POIs.append(Atom(symbol='C', position=conf_dict['liga2']['pos']))
        POIs.append(Atom(symbol='O', position=conf_dict['liga1']['N_pos']))
        POIs.append(Atom(symbol='O', position=conf_dict['liga2']['N_pos']))
        POIs.write('testing_' + str(cid) + '_POIs.xyz')
    elif type == 'ABA':
        mol.write(path='testing_' + str(cid) + '_mol.pdb', conformer=cid)
        POIs = Atoms()
        POIs.append(Atom(symbol='H', position=conf_dict['COM']))
        POIs.append(Atom(symbol='C', position=conf_dict['liga1']['pos']))
        POIs.append(Atom(symbol='C', position=conf_dict['core1']['pos']))
        POIs.append(Atom(symbol='C', position=conf_dict['liga2']['pos']))
        POIs.append(Atom(symbol='O', position=conf_dict['liga1']['N_pos']))
        POIs.append(Atom(symbol='O', position=conf_dict['liga2']['N_pos']))
        POIs.write('testing_' + str(cid) + '_POIs.xyz')


def get_binding_N_coord(molecule, conf, frag_id):
    '''Get the single (assumption) binding N in ligand building block using the
    building_block_cores of molecule.

    '''
    # assumes that the ligand block is always building block index '0'
    # get coord of possible Ns
    for i, frag in enumerate(molecule.building_block_cores(0)):
        if i == frag_id:
            frag_c = frag.GetConformer(conf)
            for atom in frag.GetAtoms():
                atom_id = atom.GetIdx()
                atom_position = frag_c.GetAtomPosition(atom_id)
                atom_position = np.array([*atom_position])
                type = frag.GetAtomWithIdx(atom_id).GetSymbol()
                if type == 'N':
                    return atom_position


def bb_center_of_mass(molecule, conf, bb_idx, frag_id):
    '''Get the center_of_mass of the core of the building block with bb_idx
    and frag_id post building.

    '''
    for i, frag in enumerate(molecule.building_block_cores(bb_idx)):
        if i == frag_id:
            center = np.array([0., 0., 0.])
            total_mass = 0.
            frag_c = frag.GetConformer(conf)
            if frag.GetAtoms():
                for atom in frag.GetAtoms():
                    atom_id = atom.GetIdx()
                    atom_position = frag_c.GetAtomPosition(atom_id)
                    atom_position = np.array([*atom_position])
                    mass = frag.GetAtomWithIdx(atom_id).GetMass()
                    total_mass += mass
                    center += mass * atom_position
                return np.divide(center, total_mass)
            else:
                return None


def get_geometrical_properties(mol, cids, type):
    '''Calculate the geometrical properties of all conformers

    '''
    # new attribute for mol
    mol.geom_prop = {}
    for cid in cids:
        # dictinary per conformer
        conf_dict = {}
        conf_dict['COM'] = mol.center_of_mass(conformer=cid)
        if type == 'ABCBA':
            # first binder -- actual ordering of BB does not matter if linker
            # information is not used
            conf_dict['liga1'] = {
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=0, frag_id=0),
                'N_pos': get_binding_N_coord(molecule=mol, conf=cid,
                                             frag_id=0)}
            # second binder -- actual ordering of BB does not matter if linker
            # information is not used
            conf_dict['liga2'] = {
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=0, frag_id=1),
                'N_pos': get_binding_N_coord(molecule=mol, conf=cid,
                                             frag_id=1)}
            # first linker -- actual ordering of BB does matter if linker
            # information is used -- currently is not
            conf_dict['link1'] = {
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=1, frag_id=0)}
            # second linker -- actual ordering of BB does matter if linker
            # information is used -- currently is not
            conf_dict['link2'] = {
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=1, frag_id=1)}
            # core
            conf_dict['core1'] = {
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=2, frag_id=0)}
        elif type == 'ABA':
            # first binder -- actual ordering of BB does not matter
            conf_dict['liga1'] = {
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=0, frag_id=0),
                'N_pos': get_binding_N_coord(molecule=mol, conf=cid,
                                             frag_id=0)}
            # second binder -- actual ordering of BB does not matter
            conf_dict['liga2'] = {
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=0, frag_id=1),
                'N_pos': get_binding_N_coord(molecule=mol, conf=cid,
                                             frag_id=1)}
            # core
            conf_dict['core1'] = {
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=1, frag_id=0)}
        mol.geom_prop[cid] = conf_dict
        # output for viz
        if False:
            # if True:
            visualize_atoms(mol, conf_dict, cid, type)

        # from the positions collected:
        # determine whether the N's are both pointing in the desired direction
        # if so,
        # determine N-N vector -> NN_v
        # determine binder core - N vectors -> BCN_1, BCN_2
        # determine binder core - binder core vector -> BCBC_v
        # calculate the angle between:
        # -- BCBC_v and BCN_i
        # -- NN_v and BCN_i (make sure the origin is not important)

        # to determine if the N's are point the right way:
        # calculate the dihedral between N1 - liga1_com - liga2_com - N2
        NBBN_dihedral = get_dihedral(pt1=mol.geom_prop[cid]['liga1']['N_pos'],
                                     pt2=mol.geom_prop[cid]['liga1']['pos'],
                                     pt3=mol.geom_prop[cid]['liga2']['pos'],
                                     pt4=mol.geom_prop[cid]['liga2']['N_pos'])
        # if the absolute value of this dihedral > some tolerance,
        # skip conformer
        mol.geom_prop[cid]['skip'] = False
        if abs(NBBN_dihedral) > 20:
            mol.geom_prop[cid]['skip'] = True
            continue
        # print(cid, NBBN_dihedral)
        mol.geom_prop[cid]['NN_v'] = mol.geom_prop[cid]['liga1']['N_pos'] \
                                     - mol.geom_prop[cid]['liga2']['N_pos']
        mol.geom_prop[cid]['BCBC_v'] = mol.geom_prop[cid]['liga1']['pos'] \
                                       - mol.geom_prop[cid]['liga2']['pos']
        mol.geom_prop[cid]['BCN_1'] = mol.geom_prop[cid]['liga1']['pos'] \
                                     - mol.geom_prop[cid]['liga1']['N_pos']
        mol.geom_prop[cid]['BCN_2'] = mol.geom_prop[cid]['liga2']['pos'] \
                                     - mol.geom_prop[cid]['liga2']['N_pos']

        # get desired angles in radian
        # negative signs applied based on the direction of vectors defined
        # above - not dependance on ordering of BB placement in stk
        mol.geom_prop[cid]['BCBC_BCN_1'] = np.degrees(
            angle_between(mol.geom_prop[cid]['BCN_1'],
                          mol.geom_prop[cid]['BCBC_v']))
        mol.geom_prop[cid]['BCBC_BCN_2'] = np.degrees(
            angle_between(mol.geom_prop[cid]['BCN_2'],
                          -mol.geom_prop[cid]['BCBC_v']))
        mol.geom_prop[cid]['NN_BCN_1'] = np.degrees(
            angle_between(mol.geom_prop[cid]['BCN_1'],
                          mol.geom_prop[cid]['NN_v']))
        mol.geom_prop[cid]['NN_BCN_2'] = np.degrees(
            angle_between(mol.geom_prop[cid]['BCN_2'],
                          -mol.geom_prop[cid]['NN_v']))

    return mol


def get_molecule(type, popns, pop_ids, N=1, mole_dir='./'):
    '''Get N conformers of a coordination cage ligand molecule using
    stk polymer function. Molecule undergoes RDKIT ETKDG conformer search and
    optimization with UFF.

    Keyword Arguments:
        type (str) - type of ligand, ABCBA or ABA
        popns (list of stk.Populations) - core, ligand and linker
            stk.Populations
        pop_ids (tuple of int) - population indices to use from core, ligand and
            linker populations
        N (int) - number of conformers to build
        mole_dir (str) - directory to save molecule to

    Returns:
        cids (list) - list of conformer IDs
        molecule (stk.Polymer) - polymer molecule with N conformers
            in molecule.mol

    '''
    core_item = popns[0][pop_ids[0]]
    liga_item = popns[1][pop_ids[1]]
    link_item = popns[2][pop_ids[2]]
    if type == 'ABCBA':
        molecule = build_ABCBA(core=core_item,
                               liga=liga_item,
                               link=link_item)
        prefix = core_item.name + '_'
        prefix += liga_item.name + '_'
        prefix += link_item.name
    elif type == 'ABA':
        molecule = build_ABA(core=core_item,
                             liga=liga_item)
        prefix = core_item.name + '_'
        prefix += liga_item.name

    # output as built
    json_file = prefix + '_' + type + '.json'
    molecule.dump(join(mole_dir, json_file))
    mol_file = prefix + '_' + type + '.mol'
    molecule.write(join(mole_dir, mol_file))
    # clean molecule
    rdkit_ETKDG(molecule)
    # output energy minimized
    json_file = prefix + '_' + type + '_opt.json'
    molecule.dump(join(mole_dir, json_file))
    mol_file = prefix + '_' + type + '_opt.mol'
    molecule.write(join(mole_dir, mol_file))
    # make N conformers of the polymer molecule
    cids = Chem.EmbedMultipleConfs(mol=molecule.mol, numConfs=N,
                                   params=Chem.ETKDG())
    # output each conformer to 3D structure if desired
    # for cid in cids:
    #     print(cid)
    #     mol_file = prefix + '_' + type + '_' + str(cid) + '_opt.mol'
    #     molecule.write(path=join(mole_dir, mol_file), conformer=cid)
    return cids, molecule


def main():
    if (not len(sys.argv) == 2):
        print("""
    Usage: molecule_builing.py N
        N (int) - number of conformers to use
        """)
        sys.exit()
    else:
        N = int(sys.argv[1])

    proj_dir = '/home/atarzia/projects/ligand_combiner/'
    core_dir = proj_dir + 'cores/'
    liga_dir = proj_dir + 'ligands/'
    link_dir = proj_dir + 'linkers/'
    mole_dir = proj_dir + 'molecules/'

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
        if inversion_flag is 't':
            i.invertable = True
        else:
            i.invertable = False

    print(core_pop)
    # this is the resultant molecule population
    molecule_pop = Population()
    for i, core in enumerate(core_pop):
        for j, liga in enumerate(liga_pop):
            for k, link in enumerate(link_pop):
                # build ABCBA molecule
                pop_ids = (i, j, k)  # core, ligand, linker
                print(pop_ids)
                ABCBA_confs, ABCBA_molecule = get_molecule(
                    type='ABCBA',
                    popns=(core_pop, liga_pop, link_pop),
                    pop_ids=pop_ids, N=N,
                    mole_dir=mole_dir)
                # get properties - save to molecule as attribute
                ABCBA_molecule = get_geometrical_properties(mol=ABCBA_molecule,
                                                            cids=ABCBA_confs,
                                                            type='ABCBA')
                molecule_pop.members.append(ABCBA_molecule)
                if ABCBA_molecule.geom_prop[0]['skip'] is False:
                    print(ABCBA_molecule.geom_prop[0]['NN_BCN_1'],
                          ABCBA_molecule.geom_prop[0]['NN_BCN_2'])
                # break
            # build ABA molecule
            ABA_confs, ABA_molecule = get_molecule(
                type='ABA',
                popns=(core_pop, liga_pop, link_pop),
                pop_ids=pop_ids, N=N,
                mole_dir=mole_dir)
            # get properties - save to molecule as attribute
            ABA_molecule = get_geometrical_properties(mol=ABA_molecule,
                                                      cids=ABA_confs,
                                                      type='ABA')
            molecule_pop.members.append(ABA_molecule)
            if ABA_molecule.geom_prop[0]['skip'] is False:
                print(ABA_molecule.geom_prop[0]['NN_BCN_1'],
                      ABA_molecule.geom_prop[0]['NN_BCN_2'])
            break
        break

    # draw 2D representation for all built molecules
    mol_list = []
    for poly in molecule_pop:
        MOL = Chem.MolFromSmiles(Chem.MolToSmiles(poly.mol))
        mol_list.append(MOL)
    mol_list2grid(molecules=mol_list, filename='built_molecules',
                  mol_per_row=3, maxrows=3, subImgSize=(200, 200))

    # find matching pairs in molecule DB


if __name__ == "__main__":
    main()
