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
from rdkit.Chem import AllChem as Chem
from os.path import join
from stk import rdkit_ETKDG
sys.path.insert(0, '/home/atarzia/thesource/')
from stk_functions import build_population, build_ABCBA, build_ABA


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
            # first binder
            conf_dict['liga1'] = {
                # frag_id is hardcoded under the assumption that
                # the ordering is maintained from Polymer buidling
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=0, frag_id=0),
                'N_pos': get_binding_N_coord(molecule=mol, conf=cid,
                                             frag_id=0)}
            # second binder
            conf_dict['liga2'] = {
                # frag_id is hardcoded under the assumption that
                # the ordering is maintained from Polymer buidling
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=0, frag_id=1),
                'N_pos': get_binding_N_coord(molecule=mol, conf=cid,
                                             frag_id=1)}
            # first linker
            conf_dict['link1'] = {
                # frag_id is hardcoded under the assumption that
                # the ordering is maintained from Polymer buidling
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=1, frag_id=0)}
            # second linker
            conf_dict['link2'] = {
                # frag_id is hardcoded under the assumption that
                # the ordering is maintained from Polymer buidling
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=1, frag_id=1)}
            # core
            conf_dict['core1'] = {
                # frag_id is hardcoded under the assumption that
                # the ordering is maintained from Polymer buidling
                # there is only one core, so this is fine.
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=2, frag_id=0)}
        elif type == 'ABA':
            # first binder
            conf_dict['liga1'] = {
                # frag_id is hardcoded under the assumption that
                # the ordering is maintained from Polymer buidling
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=0, frag_id=0),
                'N_pos': get_binding_N_coord(molecule=mol, conf=cid,
                                             frag_id=0)}
            # second binder
            conf_dict['liga2'] = {
                # frag_id is hardcoded under the assumption that
                # the ordering is maintained from Polymer buidling
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=0, frag_id=1),
                'N_pos': get_binding_N_coord(molecule=mol, conf=cid,
                                             frag_id=1)}
            # core
            conf_dict['core1'] = {
                # frag_id is hardcoded under the assumption that
                # the ordering is maintained from Polymer buidling
                # there is only one core, so this is fine.
                'pos': bb_center_of_mass(molecule=mol, conf=cid,
                                         bb_idx=1, frag_id=0)}
        mol.geom_prop[cid] = conf_dict
        print(cid, mol.geom_prop[cid])

        # output for viz
        if type == 'ABCBA':
            mol.write(path='testing_'+str(cid)+'_mol.pdb', conformer=cid)
            POIs = Atoms()
            POIs.append(Atom(symbol='H', position=conf_dict['COM']))
            POIs.append(Atom(symbol='C', position=conf_dict['liga1']['pos']))
            # POIs.append(Atom(symbol='C', position=conf_dict['link1']['pos']))
            POIs.append(Atom(symbol='C', position=conf_dict['core1']['pos']))
            # POIs.append(Atom(symbol='C', position=conf_dict['link2']['pos']))
            POIs.append(Atom(symbol='C', position=conf_dict['liga2']['pos']))
            POIs.append(Atom(symbol='O', position=conf_dict['liga1']['N_pos']))
            POIs.append(Atom(symbol='O', position=conf_dict['liga2']['N_pos']))
            POIs.write('testing_'+str(cid)+'_POIs.xyz')
        elif type == 'ABA':
            mol.write(path='testing_'+str(cid)+'_mol.pdb', conformer=cid)
            POIs = Atoms()
            POIs.append(Atom(symbol='H', position=conf_dict['COM']))
            POIs.append(Atom(symbol='C', position=conf_dict['liga1']['pos']))
            POIs.append(Atom(symbol='C', position=conf_dict['core1']['pos']))
            POIs.append(Atom(symbol='C', position=conf_dict['liga2']['pos']))
            POIs.append(Atom(symbol='O', position=conf_dict['liga1']['N_pos']))
            POIs.append(Atom(symbol='O', position=conf_dict['liga2']['N_pos']))
            POIs.write('testing_'+str(cid)+'_POIs.xyz')
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
    if (not len(sys.argv) == 1):
        print("""
    Usage: molecule_builing.py
        """)
        sys.exit()
    else:
        # CIF = sys.argv[1]
        pass

    proj_dir = '/home/atarzia/projects/ligand_combiner/'
    core_dir = proj_dir + 'cores/'
    liga_dir = proj_dir + 'ligands/'
    link_dir = proj_dir + 'linkers/'
    mole_dir = proj_dir + 'molecules/'

    # build molecule populations
    core_pop = build_population(directory=core_dir, fgs=['bromine'])
    liga_pop = build_population(directory=liga_dir, fgs=['bromine'])
    link_pop = build_population(directory=link_dir, fgs=['bromine'])

    # build ABCBA molecule
    pop_ids = (0, 0, 0)  # core, ligand, linker
    ABCBA_confs, ABCBA_molecule = get_molecule(type='ABCBA',
                                               popns=(core_pop, liga_pop, link_pop),
                                               pop_ids=pop_ids, N=5,
                                               mole_dir=mole_dir)
    # get properties - save to molecule as attribute
    ABCBA_molecule = get_geometrical_properties(mol=ABCBA_molecule,
                                                cids=ABCBA_confs,
                                                type='ABCBA')
    # build ABA molecule
    pop_ids = (0, 0, 0)  # core, ligand, linker
    ABA_confs, ABA_molecule = get_molecule(type='ABA',
                                           popns=(core_pop, liga_pop, link_pop),
                                           pop_ids=pop_ids, N=5,
                                           mole_dir=mole_dir)
    # get properties - save to molecule as attribute
    ABA_molecule = get_geometrical_properties(mol=ABA_molecule,
                                              cids=ABA_confs,
                                              type='ABA')


if __name__ == "__main__":
    main()
