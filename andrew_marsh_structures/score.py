from scscore.scscore.standalone_model_numpy import SCScorer
from scscore.utils.SA_Score.sascorer import calculateScore as SAScore
from glob import iglob
from rdkit.Chem import AllChem as rdkit
import joblib
import os
from rdkit.Chem import DataStructs
import DataStructs
import numpy as np
from collections import defaultdict

# Create a list of all the paths to mol files. Probably best to use iglob for this.
molecule_paths = [mol for mol in iglob(#Enter path here)]

def score_molecules(molecule_paths):
    """
    Args: list - str: paths to mol files.
    Returns: dictionary - key: str:Mol in inchi. value: float/int: Scores for each mol file.
    """
    scores = defaultdict(list)
    scs_model = SCScorer()
    scs_model.restore('scscore/models/full_reaxys_model_1024bool/model.ckpt-10654.as_numpy.json.gz')
    rf_model = joblib.load('rf_model.joblib')
    for idx, mol in enumerate(molecule_paths):
        smi = rdkit.MolToSmiles(rdkit.MolFromMolFile(mol))
        (smi, sco) = scs_model.get_score_from_smi(smi)
        # SA_score prediction
        SA_score = (SAScore(rdkit.MolFromSmiles(smi)))
        #SCS_score prediction
        SCS_score = sco
        # Converts the mol file to rdkit mol, then to fingerprint (vector length 1x1024) and makes the prediction using random forest model.
        ML_model = rf_model.predict([mol2fp(rdkit.MolFromMolFile(mol))])
        scores[smi].append((SCS_score, SA_score, ML_model[0]))
    return scores

# In this dictionary the keys are smiles for the molecule athe values are the SA_Score the SCS_Score and ML_model preduictions respectively. The ML model outputs a 1 if synthesisable and a 0 if it isn't.
score = score_molecules(molecule_paths)
score

# Convert rdkit mol to fingerprint
def mol2fp(mol,nBits=1024):
    fp = rdkit.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    return fp
