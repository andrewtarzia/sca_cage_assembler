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


def get_atom_coord(molecule, conf, atom_type):
    '''Get the cartesian coordinates of 1 atom of type 'atom_type' in confomer
    'conf' in 'molecule'.

    '''
    # Get the conformer from the rdkit instance.
    conformer = molecule.GetConformer(conf)
    # check that there is only one atom of type atom_type
    count_type = 0
    for atom in molecule.GetAtoms():
        if atom.GetSymbol() == atom_type:
            des_atom_idx = atom.GetIdx()
            count_type += 1
    if count_type > 1:
        raise Exception('more then one {} in molecule'.format(atom_type))

    # get coordinates
    coord = conformer.GetAtomPosition(des_atom_idx)
    return np.array([*coord])


def get_binding_N_coord(molecule, conf, bb):
    '''Get the single (assumption) binding N in a building block bb using the
    building_block_cores of molecule.

    '''
    # assumes that the ligand block is always building block index '0'
    # get coord of possible Ns
    N_coords = [get_atom_coord(molecule=i, conf=conf, atom_type='N')
                for i in molecule.building_block_cores(0)]
    print('n coords', N_coords)
    # output the coord closest to bb COM
    dist_to_COM = []
    cent = bb.center_of_mass()
    print('COM', cent)
    for coord in N_coords:
        print('c', coord)
        # calculate catesian distance
        distance = np.linalg.norm(coord - cent)
        print('d', distance)
        dist_to_COM.append(distance)
    closest_N_coord = N_coords[dist_to_COM.index(min(dist_to_COM))]
    print('closest n', closest_N_coord)
    return closest_N_coord


def get_geometrical_properties(mol, type, cids):
    '''Calculate the geometrical properties of all conformers

    '''
    # new attribute for mol
    mol.geom_prop = {}
    for cid in cids:

        print([i.center_of_mass(conformer=cid) for i in mol.building_blocks])


        # sys.exit()

        # dictinary per conformer
        conf_dict = {}
        # type dependant dictionary - check the number of BB is correct
        if type == 'ABCBA' and len(mol.building_blocks) == 5:
            # add checks that file name matches the expected ordering
            if 'lig_' in mol.building_blocks[0].file:
                print('lig1')
                conf_dict['liga1'] = {
                    'pos': mol.building_blocks[0].center_of_mass(conformer=cid),
                    'N_pos': get_binding_N_coord(molecule=mol, conf=cid,
                                                 bb=mol.building_blocks[0])}
            else:
                raise Exception('building block ordering is off.')
            if 'link_' in mol.building_blocks[1].file:
                conf_dict['link1'] = {
                    'pos': mol.building_blocks[1].center_of_mass(conformer=cid)}
            else:
                raise Exception('building block ordering is off.')
            if 'core_' in mol.building_blocks[2].file:
                conf_dict['core1'] = {
                    'pos': mol.building_blocks[2].center_of_mass(conformer=cid)}
            else:
                raise Exception('building block ordering is off.')
            if 'link_' in mol.building_blocks[3].file:
                conf_dict['link2'] = {
                    'pos': mol.building_blocks[3].center_of_mass(conformer=cid)}
            else:
                raise Exception('building block ordering is off.')
            if 'lig_' in mol.building_blocks[4].file:
                print('lig2')
                conf_dict['liga2'] = {
                    'pos': mol.building_blocks[4].center_of_mass(conformer=cid),
                    'N_pos': get_binding_N_coord(molecule=mol, conf=cid,
                                                 bb=mol.building_blocks[4])}
            else:
                raise Exception('building block ordering is off.')
        elif type == 'ABA' and len(mol.building_blocks) == 3:
            if 'lig_' in mol.building_blocks[0].file:
                print('lig1')
                conf_dict['liga1'] = {
                    'pos': mol.building_blocks[0].center_of_mass(conformer=cid),
                    'N_pos': get_binding_N_coord(molecule=mol, conf=cid,
                                                 bb=mol.building_blocks[0])}
            else:
                raise Exception('building block ordering is off.')
            if 'core_' in mol.building_blocks[1].file:
                conf_dict['core1'] = {
                    'pos': mol.building_blocks[1].center_of_mass(conformer=cid)}
            else:
                raise Exception('building block ordering is off.')
            if 'lig_' in mol.building_blocks[2].file:
                print('lig2')
                conf_dict['liga2'] = {
                    'pos': mol.building_blocks[2].center_of_mass(conformer=cid),
                    'N_pos': get_binding_N_coord(molecule=mol, conf=cid,
                                                 bb=mol.building_blocks[2])}
            else:
                raise Exception('building block ordering is off.')
        else:
            raise Exception('something went wrong here, either type or build process')
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
                                                type='ABCBA',
                                                cids=ABCBA_confs)
    # build ABA molecule
    pop_ids = (0, 0, 0)  # core, ligand, linker
    ABA_confs, ABA_molecule = get_molecule(type='ABA',
                                           popns=(core_pop, liga_pop, link_pop),
                                           pop_ids=pop_ids, N=5,
                                           mole_dir=mole_dir)
    # get properties - save to molecule as attribute
    ABA_molecule = get_geometrical_properties(mol=ABA_molecule,
                                              type='ABA',
                                              cids=ABA_confs)


if __name__ == "__main__":
    main()
