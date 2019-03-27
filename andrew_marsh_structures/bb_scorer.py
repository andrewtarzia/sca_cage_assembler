#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to obtain synthetic accessibility of a DB of mol files using 3 methods.

Author: Andrew Tarzia

Code borrowed from score.py by Steven Bennett.

Date Created: 27 Feb 2019
"""

import sys
import glob
from scscore.scscore.standalone_model_numpy import SCScorer
from scscore.utils.SA_Score.sascorer import calculateScore as SAScore
from rdkit.Chem import AllChem as Chem
import joblib
import os
import numpy as np


def mol2fp(mol,nBits=1024):
    '''Get Morgan fingerprint of RDKIT Molecule.

    '''
    fp = rdkit.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    return fp


def write_output(file, name, SAscore, SCscore, ML):
    '''Write molecule properties to output file.

    '''
    with open(file, 'a') as f:
        f.write(name+','+str(SAscore)+','+str(SCscore)+','+str(ML)+'\n')


def get_SA(smi):
    '''Get SAscore for a given smiles.

    '''
    SA = (SAScore(rdkit.MolFromSmiles(smi)))
    return SA


def get_SC(smi, model):
    '''Get SAscore for a given smiles.

    '''
    _, SC = model.get_score_from_smi(smi)
    return SC


def get_ML(smi, model):
    '''Get SAscore for a given smiles.

    '''
    ML = model.predict([mol2fp(rdkit.MolFromSmiles(smi))])
    return ML[0]


def main():
    """Run script.

    """
    if (not len(sys.argv) == 2):
        print("""
Usage: bb_scorer.py output
    output: file to output results""")
        sys.exit()
    else:
        output = sys.argv[1]
    src_dir = '/home/atarzia/projects/andrew_marsh_structures/src/'

    # initiliaze Lukas ML model
    ML_model = joblib.load(src_dir+'rf_model.joblib')
    # initialize SCscore ML model
    SC_model = SCScorer()
    SC_model.restore(src_dir+'scscore/models/full_reaxys_model_1024bool/model.ckpt-10654.as_numpy.json.gz')
    # init output file
    with open(output, 'w') as f:
        f.write('name,SA,SC,ML\n')
    mol_list = glob.glob('*.mol')
    print('doing', len(mol_list), 'molecules')
    for mol in mol_list:
        name = mol.replace('.mol', '')
        # read into rdkit -> SMILES
        smi = Chem.MolToSmiles(Chem.MolFromMolFile(mol))
        # get SA score (ertl)
        SAscore = get_SA(smi=smi)
        # get SC score ()
        SCscore = get_SC(smi=smi, model=SC_model)
        # get ML score (Lukas)
        ML = get_ML(smi=smi, model=ML_model)
        # output
        write_output(file=output, name=name, SAscore=SAscore,
                     SCscore=SCscore, ML=ML)
    print('done')


if __name__ == "__main__":
    main()
