#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for containing functions for ligand_combiner

Author: Andrew Tarzia

Date Created: 17 Apr 2019

"""

import sys
from ase.atoms import Atoms, Atom
import numpy as np
from rdkit.Chem import AllChem as Chem
from os.path import join
from stk import rdkit_ETKDG
sys.path.insert(0, '/home/atarzia/thesource/')
from stk_functions import build_ABCBA, build_ABA
from calculations import get_dihedral, angle_between


class Combination:
    '''Class defining a pair of linkers and their geometrical properties

    '''
    def __init__(self, mol1, mol2, conf1, conf2):
        self.mol1 = mol1
        self.mol2 = mol2
        self.conf1 = conf1
        self.conf2 = conf2

    def get_angle_deviations(self):
        '''Get the closeness of test angles to 180 degrees

        '''
        self.angle1_deviation = abs(180 - (self.p1_angle1 + self.p2_angle1))
        self.angle2_deviation = abs(180 - (self.p1_angle2 + self.p2_angle2))

    def get_N_Pd_lengths_deviation(self):
        '''Get the difference of the LHS and RHS that are defined by the
        N-Pd-N distance

        '''
        self.NPdN_difference = abs(self.extender_V_LHS-self.extender_V_RHS)

    def test_N_N_lengths(self):
        '''Test that the larger linker has a longer NN distance than the shorter
        linker.

        '''
        if self.NN_dist1 > self.NN_dist2:
            True
        else:
            return False


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
                return np.array([0., 0., 0.])
    return np.array([0., 0., 0.])


def get_geometrical_properties(mol, cids, type):
    '''Calculate the geometrical properties of all conformers

    '''
    # new attribute for mol
    mol.geom_prop = {}
    for cid in cids:
        # dictinary per conformer
        conf_dict = {}
        conf_dict['COM'] = mol.center_of_mass(conformer=cid).tolist()
        if type == 'ABCBA':
            # first binder -- actual ordering of BB does not matter if linker
            # information is not used
            conf_dict['liga1'] = {
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=0, frag_id=0).tolist(),
                'N_pos': get_binding_N_coord(molecule=mol, conf=cid,
                                             frag_id=0).tolist(),
                'CC_pos': get_binding_N_CC_coord(molecule=mol, conf=cid,
                                                 frag_id=0).tolist()}
            # second binder -- actual ordering of BB does not matter if linker
            # information is not used
            conf_dict['liga2'] = {
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=0, frag_id=1).tolist(),
                'N_pos': get_binding_N_coord(molecule=mol, conf=cid,
                                             frag_id=1).tolist(),
                'CC_pos': get_binding_N_CC_coord(molecule=mol, conf=cid,
                                                 frag_id=1).tolist()}
            # first linker -- actual ordering of BB does matter if linker
            # information is used -- currently is not
            conf_dict['link1'] = {
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=1, frag_id=0).tolist()}
            # second linker -- actual ordering of BB does matter if linker
            # information is used -- currently is not
            conf_dict['link2'] = {
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=1, frag_id=1).tolist()}
            # core
            conf_dict['core1'] = {
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=2, frag_id=0).tolist()}
        elif type == 'ABA':
            # first binder -- actual ordering of BB does not matter
            conf_dict['liga1'] = {
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=0, frag_id=0).tolist(),
                'N_pos': get_binding_N_coord(molecule=mol, conf=cid,
                                             frag_id=0).tolist(),
                'CC_pos': get_binding_N_CC_coord(molecule=mol, conf=cid,
                                                 frag_id=0).tolist()}
            # second binder -- actual ordering of BB does not matter
            conf_dict['liga2'] = {
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=0, frag_id=1).tolist(),
                'N_pos': get_binding_N_coord(molecule=mol, conf=cid,
                                             frag_id=1).tolist(),
                'CC_pos': get_binding_N_CC_coord(molecule=mol, conf=cid,
                                                 frag_id=1).tolist()}
            # core
            conf_dict['core1'] = {
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=1, frag_id=0).tolist()}
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
        dihedral_lim = 20
        if abs(NBBN_dihedral) > dihedral_lim:
            mol.geom_prop[cid]['skip'] = True
            # print('skipping, NBBN dihedral > {}'.format(dihedral_lim))
            continue
        # also check that the N atoms are pointing away from the core
        # by making sure COM_core-N_pos > COM_core-liga_pos
        core_2_COM = np.asarray(mol.geom_prop[cid]['COM']) \
                     - np.asarray(mol.geom_prop[cid]['core1']['pos'])
        core_2_N1 = np.asarray(mol.geom_prop[cid]['liga1']['N_pos']) \
                     - np.asarray(mol.geom_prop[cid]['core1']['pos'])
        core_2_liga1 = np.asarray(mol.geom_prop[cid]['liga1']['pos']) \
                     - np.asarray(mol.geom_prop[cid]['core1']['pos'])
        core_2_N2 = np.asarray(mol.geom_prop[cid]['liga2']['N_pos']) \
                     - np.asarray(mol.geom_prop[cid]['core1']['pos'])
        core_2_liga2 = np.asarray(mol.geom_prop[cid]['liga2']['pos']) \
                     - np.asarray(mol.geom_prop[cid]['core1']['pos'])

        N1_core_COM_angle = angle_between(core_2_N1, core_2_COM)
        L1_core_COM_angle = angle_between(core_2_liga1, core_2_COM)
        N2_core_COM_angle = angle_between(core_2_N2, core_2_COM)
        L2_core_COM_angle = angle_between(core_2_liga2, core_2_COM)

        if N1_core_COM_angle > L1_core_COM_angle:
            if N2_core_COM_angle > L2_core_COM_angle:
                mol.geom_prop[cid]['skip'] = True
                # print("skipping, both N's point backwards")
                continue

        NN_v = np.asarray(mol.geom_prop[cid]['liga1']['N_pos']) \
                          - np.asarray(mol.geom_prop[cid]['liga2']['N_pos'])
        BCN_1 = np.asarray(mol.geom_prop[cid]['liga1']['CC_pos']) \
                           - np.asarray(mol.geom_prop[cid]['liga1']['N_pos']).tolist()
        BCN_2 = np.asarray(mol.geom_prop[cid]['liga2']['CC_pos']) \
                           - np.asarray(mol.geom_prop[cid]['liga2']['N_pos']).tolist()
        mol.geom_prop[cid]['NN_v'] = NN_v.tolist()
        mol.geom_prop[cid]['BCN_1'] = BCN_1.tolist()
        mol.geom_prop[cid]['BCN_2'] = BCN_2.tolist()
        print('{}: confomer {} passed'.format(mol.name, cid))
        # output for viz
        # if False:
        if True:
            visualize_atoms(mol, conf_dict, cid, type,
                            filename=mol.name + '_' + str(cid))

        # get desired angles in radian
        # negative signs applied based on the direction of vectors defined
        # above - not dependance on ordering of BB placement in stk
        mol.geom_prop[cid]['NN_BCN_1'] = float(np.degrees(
            angle_between(np.asarray(mol.geom_prop[cid]['BCN_1']),
                          np.asarray(mol.geom_prop[cid]['NN_v']))))
        mol.geom_prop[cid]['NN_BCN_2'] = float(np.degrees(
            angle_between(np.asarray(mol.geom_prop[cid]['BCN_2']),
                          -np.asarray(mol.geom_prop[cid]['NN_v']))))

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
    etkdg = Chem.ETKDG()
    etkdg.randomSeed = 1000
    cids = Chem.EmbedMultipleConfs(mol=molecule.mol, numConfs=N,
                                   params=etkdg)
    # output each conformer to 3D structure if desired
    # for cid in cids:
    #     print(cid)
    #     mol_file = prefix + '_' + type + '_' + str(cid) + '_opt.mol'
    #     molecule.write(path=join(mole_dir, mol_file), conformer=cid)
    return cids, molecule
