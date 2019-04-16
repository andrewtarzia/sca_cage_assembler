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
import stk
from stk import rdkit_ETKDG, Population
sys.path.insert(0, '/home/atarzia/thesource/')
from stk_functions import build_population, build_ABCBA, build_ABA
from calculations import get_dihedral, angle_between
from rdkit_functions import mol_list2grid


def atoms_2_vect(ASE, p1, p2):
    '''Append to ASE.Atoms() the interpolation between two points.

    '''
    pts = [np.linspace(p1[i], p2[i]) for i in np.arange(len(p1))]
    for i, j, k in zip(*pts):
        ASE.append(Atom(symbol='P', position=[i, j, k]))
    return ASE


def visualize_atoms(mol, conf_dict, cid, type, filename):
    '''Output files with POIs and conformer atom positions for visulization.

    '''

    mol.write(path=filename + '.pdb', conformer=cid)
    POIs = Atoms()
    POIs.append(Atom(symbol='H', position=conf_dict['COM']))
    POIs.append(Atom(symbol='C', position=conf_dict['liga1']['pos']))
    POIs.append(Atom(symbol='C', position=conf_dict['core1']['pos']))
    POIs.append(Atom(symbol='C', position=conf_dict['liga2']['pos']))
    POIs.append(Atom(symbol='Be', position=conf_dict['liga1']['CC_pos']))
    POIs.append(Atom(symbol='Be', position=conf_dict['liga2']['CC_pos']))
    POIs.append(Atom(symbol='O', position=conf_dict['liga1']['N_pos']))
    POIs.append(Atom(symbol='O', position=conf_dict['liga2']['N_pos']))
    # plot vectors as P atom
    POIs = atoms_2_vect(ASE=POIs, p1=conf_dict['liga1']['CC_pos'],
                        p2=conf_dict['liga1']['N_pos'])
    POIs = atoms_2_vect(ASE=POIs, p1=conf_dict['liga2']['CC_pos'],
                        p2=conf_dict['liga2']['N_pos'])
    POIs = atoms_2_vect(ASE=POIs, p1=conf_dict['liga1']['N_pos'],
                        p2=conf_dict['liga2']['N_pos'])
    # if type == 'ABCBA':
    #    POIs.append(Atom(symbol='C', position=conf_dict['link1']['pos']))
    #    POIs.append(Atom(symbol='C', position=conf_dict['link2']['pos']))
    POIs.write(filename + '_POIs.xyz')


def get_binding_N_CC_coord(molecule, conf, frag_id):
    '''Get the midpoint of the vector connecting the single (assumption)
    binding N in ligand building block using the building_block_cores of
    molecule.

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
                    # get neighbours of this N
                    C_ids = []
                    for atom2 in atom.GetNeighbors():
                        atom2_id = atom2.GetIdx()
                        # if they are carbons
                        if frag.GetAtomWithIdx(atom2_id).GetSymbol():
                            C_ids.append(atom2_id)
                    CC1_id = C_ids[0]
                    CC2_id = C_ids[1]
                    # get their atom positions
                    CC1_pos = frag_c.GetAtomPosition(CC1_id)
                    CC1 = np.array([*CC1_pos])
                    CC2_pos = frag_c.GetAtomPosition(CC2_id)
                    CC2 = np.array([*CC2_pos])
                    # get the vector between C neighbours
                    CC_v = CC1 - CC2
                    # get the midpoint of the vector
                    CC_midpoint = CC1 - CC_v / 2
                    return CC_midpoint


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
                                             frag_id=0),
                'CC_pos': get_binding_N_CC_coord(molecule=mol, conf=cid,
                                                 frag_id=0)}
            # second binder -- actual ordering of BB does not matter if linker
            # information is not used
            conf_dict['liga2'] = {
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=0, frag_id=1),
                'N_pos': get_binding_N_coord(molecule=mol, conf=cid,
                                             frag_id=1),
                'CC_pos': get_binding_N_CC_coord(molecule=mol, conf=cid,
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
                                             frag_id=0),
                'CC_pos': get_binding_N_CC_coord(molecule=mol, conf=cid,
                                                 frag_id=0)}
            # second binder -- actual ordering of BB does not matter
            conf_dict['liga2'] = {
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=0, frag_id=1),
                'N_pos': get_binding_N_coord(molecule=mol, conf=cid,
                                             frag_id=1),
                'CC_pos': get_binding_N_CC_coord(molecule=mol, conf=cid,
                                                 frag_id=1)}
            # core
            conf_dict['core1'] = {
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=1, frag_id=0)}
        mol.geom_prop[cid] = conf_dict
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
        mol.geom_prop[cid]['BCN_1'] = mol.geom_prop[cid]['liga1']['CC_pos'] \
                                     - mol.geom_prop[cid]['liga1']['N_pos']
        mol.geom_prop[cid]['BCN_2'] = mol.geom_prop[cid]['liga2']['CC_pos'] \
                                     - mol.geom_prop[cid]['liga2']['N_pos']

        # output for viz
        # if False:
        if True:
            visualize_atoms(mol, conf_dict, cid, type,
                            filename=mol.name + '_' + str(cid))

        # get desired angles in radian
        # negative signs applied based on the direction of vectors defined
        # above - not dependance on ordering of BB placement in stk
        mol.geom_prop[cid]['NN_BCN_1'] = np.degrees(
            angle_between(mol.geom_prop[cid]['BCN_1'],
                          mol.geom_prop[cid]['NN_v']))
        mol.geom_prop[cid]['NN_BCN_2'] = np.degrees(
            angle_between(mol.geom_prop[cid]['BCN_2'],
                          -mol.geom_prop[cid]['NN_v']))

    return mol


def get_molecule(type, popns, pop_ids, inverted, N=1, mole_dir='./'):
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
                               link=link_item,
                               flippedlink=inverted)
        prefix = core_item.name + '_'
        prefix += liga_item.name + '_'
        if inverted:
            prefix += link_item.name + 'i'
        else:
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
            for j, liga in enumerate(liga_pop):
                for k, link in enumerate(link_pop):
                    # build ABCBA molecule
                    pop_ids = (i, j, k)  # core, ligand, linker
                    print(pop_ids)
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
        core_pop.dump(join(core_dir, 'core.pop'))
        liga_pop.dump(join(liga_dir, 'ligands.pop'))
        link_pop.dump(join(link_dir, 'linkers.pop'))
        molecule_pop.dump(join(mole_dir, 'molecules.pop'))

        # draw 2D representation for all built molecules
        mol_list = []
        for poly in molecule_pop:
            MOL = Chem.MolFromSmiles(Chem.MolToSmiles(poly.mol))
            mol_list.append(MOL)
        mol_list2grid(molecules=mol_list, filename='built_molecules',
                      mol_per_row=3, maxrows=3, subImgSize=(200, 200))

    elif rebuild == 'f':
        # load in populations
        core_pop = stk.Population()
        liga_pop = stk.Population()
        link_pop = stk.Population()
        molecule_pop = stk.Population()
        core_pop = core_pop.load(path=join(core_dir, 'core.pop'),
                                 member_init=stk.Molecule.from_dict)
        liga_pop = liga_pop.load(path=join(liga_dir, 'ligands.pop'),
                                 member_init=stk.Molecule.from_dict)
        link_pop = link_pop.load(path=join(link_dir, 'linkers.pop'),
                                 member_init=stk.Molecule.from_dict)
        molecule_pop = molecule_pop.load(path=join(mole_dir, 'molecules.pop'),
                                         member_init=stk.Molecule.from_dict)

    print(molecule_pop)

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
            PROP1 = poly1.geom_prop[conf1]
            # skip conformer if dihedral meant the N's were not on the right side
            if PROP1['skip'] is True:
                continue
            for j, poly2 in enumerate(molecule_pop):
                # make sure poly1 != poly2
                if i != j:
                    for conf2 in poly2.geom_prop:
                        print(poly1.name, poly2.name)
                        PROP2 = poly2.geom_prop[conf2]
                        # skip conformer if desired
                        if PROP2['skip'] is True:
                            continue
                        # check N-N distance of poly1-conf > poly2-conf
                        NN_dist1 = np.linalg.norm(PROP1['NN_v'])
                        NN_dist2 = np.linalg.norm(PROP2['NN_v'])
                        print(NN_dist1, NN_dist2)
                        if NN_dist1 > NN_dist2:
                            # determine angles made by NN_v and NN-BC_v
                            # check that the pairs sum to 180
                            p1_angle1 = PROP1['NN_BCN_1']
                            p1_angle2 = PROP1['NN_BCN_2']
                            p2_angle1 = PROP2['NN_BCN_1']
                            p2_angle2 = PROP2['NN_BCN_2']
                            print(i, j, p1_angle1, p1_angle2, p2_angle1, p2_angle2)
                            if np.isclose(p1_angle1 + p2_angle1, 180, rtol=0, atol=angle_tol):
                                if np.isclose(p1_angle2 + p2_angle2, 180, rtol=0, atol=angle_tol):
                                    # now check that the length of the long
                                    # vector and the short vector are commensurate
                                    # with an ideal trapezoid with the given angles
                                    # i.e. the extender vector determined by the
                                    # difference of the two NN_v (LHS) matches what is
                                    # expected by trig (RHS)
                                    extender_V_LHS = (NN_dist1 - NN_dist2) / 2
                                    print(extender_V_LHS)
                                    test_angle = np.radians(180 - p2_angle1)
                                    print(test_angle)
                                    extender_V_RHS = vector_length * np.cos(test_angle)
                                    print(extender_V_RHS)
                                    tol = vector_std * np.cos(test_angle)
                                    print(tol)
                                    if np.isclose(extender_V_LHS, extender_V_RHS,
                                                  rtol=0, atol=tol):
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


if __name__ == "__main__":
    main()
