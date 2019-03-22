#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions that are useful for rdkit usage

Author: Andrew Tarzia

Date Created: 15 Mar 2019
"""

import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors, Draw, PyMol
from rdkit.Chem.Descriptors3D import NPR1, NPR2, PMI1, PMI2, PMI3
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from rdkit.Geometry import rdGeometry
from rdkit import Geometry


def calculate_all_MW(molecules):
    """Calculate the molecular weight of all molecules in DB dictionary.

    {name: SMILES}

    """
    for m, smile in molecules.items():
        # Read SMILES and add Hs
        mol = Chem.AddHs(Chem.MolFromSmiles(smile))
        MW = Descriptors.MolWt(mol)
        print(m, '---', smile, '---', 'MW =', MW, 'g/mol')


def draw_smiles_to_svg(smiles, filename):
    """Draw a single molecule to an SVG file with transparent BG.

    """
    mol = Chem.MolFromSmiles(smiles)
    # change BG to transperent
    # (https://sourceforge.net/p/rdkit/mailman/message/31637105/)
    o = DrawingOptions()
    o.bgColor = None
    Chem.Compute2DCoords(mol)
    Draw.MolToFile(mol, filename, fitImage=True, imageType='svg',
                   options=o)


def mol_list2grid(mol_dict, filename, mol_per_row, subImgSize=(200, 200)):
    '''Produce a grid of molecules in mol_list.

    mol_dict (dict) - {name: SMILES}

    '''
    mols = []
    names = []
    for name in mol_dict:
        names.append(name)
        mols.append(Chem.MolFromSmiles(mol_dict[name]))
    img = Draw.MolsToGridImage(mols, molsPerRow=mol_per_row,
                               subImgSize=subImgSize,
                               legends=names,
                               useSVG=False)
    img.save(filename+'.png')


def read_mol_txt_file(filename):
    """Function to read molecule SMILES and information from txt file.

    """
    data = pd.read_table(filename, delimiter=':')
    molecules = {}
    diameters = {}
    for i, row in data.iterrows():
        # try:
        #     name, smile, radius = line.rstrip().split(':')
        # except ValueError:
        #     print(line, 'had : in there twice, fix this naming or SMILE')
        #     print('skipped')
        name = row['molecule']
        smile = row['smile']
        diameter = row['diameter']
        molecules[name] = smile
        diameters[name] = diameter
    return data, molecules, diameters


def get_inertial_prop(mol, cids):
    """Get inertial 3D descriptors for all conformers in mol.

    """
    # ratio 1 is I1/I3
    # ratio 2 is I2/I3
    sml_PMI, mid_PMI, lge_PMI = [], [], []
    ratio_1_, ratio_2_ = [], []
    for cid in cids:
        sml_PMI.append(PMI1(mol, confId=cid))
        mid_PMI.append(PMI2(mol, confId=cid))
        lge_PMI.append(PMI3(mol, confId=cid))
        ratio_1_.append(NPR1(mol, confId=cid))
        ratio_2_.append(NPR2(mol, confId=cid))

    return sml_PMI, mid_PMI, lge_PMI, ratio_1_, ratio_2_


def get_COMs(mol, cids):
    """Get COM of all conformers of mol.

    Code from:
    https://iwatobipen.wordpress.com/2016/08/16/scoring-3d-diversity-using-rdkit-rdkit/

    """
    coms = []
    numatoms = mol.GetNumAtoms()
    for confId in range(len(cids)):
        # print('conf:', confId)
        # print('number of atoms:', numatoms)
        conf = mol.GetConformer(confId)
        coords = np.array([list(conf.GetAtomPosition(atmidx)) for atmidx in range(numatoms)])
        # print('coords:')
        # print(coords)
        atoms = [atom for atom in mol.GetAtoms()]
        mass = Descriptors.MolWt(mol)
        # print('mass:', mass)
        centre_of_mass = np.array(np.sum(atoms[i].GetMass() * coords[i] for i in range(numatoms))) / mass
        # print(centre_of_mass)
        coms.append(centre_of_mass)

    return coms


def def_point(x, y, z):
    """Define a 3D point in RDKIT

    """
    point = rdGeometry.Point3D()
    point.x = x
    point.y = y
    point.z = z

    return point


def smiles2conformers(smiles, N=10, optimize=True):
    '''Convert smiles string to N conformers.

    Keyword Arguments:
        smiles (str) - smiles string for molecule
        N (int) - number of conformers to generate using the ETKDG algorithm
        optimize (bool) - flag for UFF optimization (default=True)

    Returns:
        mol (RDKit molecule ::class::) - contains N conformers
    '''
    # Read SMILES and add Hs
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print('RDKit error for', smiles)
        return None
    mol = Chem.AddHs(mol)
    # try based on RuntimeError from RDKit
    try:
        # 2D to 3D with multiple conformers
        cids = Chem.EmbedMultipleConfs(mol=mol, numConfs=N,
                                       useExpTorsionAnglePrefs=True,
                                       useBasicKnowledge=True)
        # quick UFF optimize
        for cid in cids:
            Chem.UFFOptimizeMolecule(mol, confId=cid)
    except RuntimeError:
        print('RDKit error for', smiles)
        return None
    return mol
