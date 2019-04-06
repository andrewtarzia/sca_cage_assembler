#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build all molecules from DB of core, ligands and linkers.

Author: Andrew Tarzia

Date Created: 04 Apr 2019

"""

import sys
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
def get_geometrical_properties(mol, type, cids):
    '''Calculate the geometrical properties of all conformers

    '''
    for cid in cids:
        print(cid)
        targ_mol = mol.mol.GetConformer(cid)
        bb = mol.building_blocks[0]
        print([i.func_group_infos[0].name for i in mol.building_blocks])
        print(mol.fg_centroid(fg_id=1, conformer=cid))
        print(mol.bonder_ids)
        print(targ_mol)
        sys.exit()


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
