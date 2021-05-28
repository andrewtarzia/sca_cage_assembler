#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to extract cage molecule from CIF.

Author: Andrew Tarzia

Date Created: 15 Mar 2019
"""
import sys
from rdkit.Chem import AllChem as Chem


def get_guests(paper):
    guests = {
        'smulders2013':
        {
            'acetone': 'CC(=O)C',
            'pyridine': 'C1=CC=NC=C1',
            'pyridazine': 'C1=CC=NN=C1',
            'pyrimidine': 'C1=CN=CN=C1',
            'pyrazine': 'C1=CN=CC=N1',
            'tetrahydropyran': 'C1CCOCC1',
            '1,3-dioxane': 'C1COCOC1',
            '1,4-dioxane': 'C1COCCO1',
            '1,3,5-trioxane': 'C1OCOCO1',
            'tetrahydrofuran': 'C1CCOC1',
            'furan': 'C1=COC=C1',
            'CH2Cl2': 'C(Cl)Cl',
            'CHCl3': 'C(Cl)(Cl)Cl',
            'CCl4': 'C(Cl)(Cl)(Cl)Cl',
            'SF6': 'FS(F)(F)(F)(F)F',
            'methylcyclopentane': 'CC1CCCC1',
            'cyclohexane': 'C1CCCCC1',
            'methylcyclohexane': 'CC1CCCCC1',
            'cycloheptane': 'C1CCCCCC1',
            'cyclohexene': 'C1CCC=CC1',
            '1,3-cyclohexadiene': 'C1CC=CC=C1',
            'benzene': 'C1=CC=CC=C1',
            'fluorobenzene': 'C1=CC=C(C=C1)F',
            'toluene': 'CC1=CC=CC=C1',
        },
        'bolliger2014':
        {
            'acetonitrile': 'CC#N',
            'water': 'O',
            '2-formylpyridine': 'C1=CC=NC(=C1)C=O',
            'diethylether': 'CCOCC',
            'pyridine': 'C1=CC=NC=C1',
            'tetrahydrofuran': 'C1CCOC1',
            'mesitylene': 'CC1=CC(=CC(=C1)C)C',
            'cyclopentanone': 'C1CCC(=O)C1',
            'm-xylene': 'CC1=CC(=CC=C1)C',
            'camphor': 'CC1(C2CCC1(C(=O)C2)C)C',
            'hexafluorobenzene': 'C1(=C(C(=C(C(=C1F)F)F)F)F)F',
            'cyclopentane': 'C1CCCC1',
            '1-adamantylmethanol': 'C1C2CC3CC1CC(C2)(C3)CO',
            '2-methylthiophene': 'CC1=CC=CS1',
            'cyclohexane': 'C1CCCCC1',
            'hexane': 'CCCCCC',
            '1,5-cyclooctadiene': "C1/C=C\CC/C=C\C1",
            '2-hexylthiophene': 'CCCCCCC1=CC=CS1',
            'trans-decalin': '	[H][C@@]12CCCC[C@@]1([H])CCCC2',
            'cis-decalin': '[H][C@@]12CCCC[C@]1([H])CCCC2',
            'cyclooctane': 'C1CCCCCCC1',
            'adamantane': 'C1C2CC3CC1CC(C2)C3',
            'R-limonene': 'CC1=CC[C@@H](CC1)C(=C)C',
            'naphthalene': 'C1=CC=C2C=CC=CC2=C1',
            '1,3,5-Triethylbenzene': 'CCC1=CC(=CC(=C1)CC)CC',
            '8-phenyloctanol': 'C1=CC=C(C=C1)CCCCCCCCO',
            'tetraphenylmethane': 'C1=CC=C(C=C1)C(C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4',
            'diphenylether': 'C1=CC=C(C=C1)OC2=CC=CC=C2',
            'biphenyl': 'C1=CC=C(C=C1)C2=CC=CC=C2',
            'phenanthrene': 'C1=CC=C2C(=C1)C=CC3=CC=CC=C32',
            'pyrene': 'C1=CC2=C3C(=C1)C=CC4=CC=CC(=C43)C=C2',
            'fluorene': 'C1C2=CC=CC=C2C3=CC=CC=C31',
            'anthracene': 'C1=CC=C2C=C3C=CC=CC3=CC2=C1',
            'D-glucose': 'C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O',
            'cyclododecane': 'C1CCCCCCCCCCC1',
            '18-crown-6': 'C1COCCOCCOCCOCCOCCO1',
            '15-crown-5': 'C1COCCOCCOCCOCCO1',
            'D-leucine': 'CC(C)C[C@H](C(=O)O)N',
        }
    }
    if paper not in guests:
        print('paper must be smulders2013 or bolliger2014')
        sys.exit('exitting.')
    return guests[paper]


if __name__ == "__main__":
    if (not len(sys.argv) == 2):
        print("""
Usage: guests.py paper
    paper: paper whose guests you want (smulders2013 or bolliger2014)
    """)
        sys.exit()
    else:
        paper = sys.argv[1]

    guests = get_guests(paper=paper)
    print(guests)
    for guest in guests:
        if guest in ['SF6']:
            # problematic for RDKit, so skip them
            # built manually
            continue
        # get RDKIT mol -- get a 3D conformer using ETKDG --
        # energy minimize with UFF
        mol = smiles2conformers(smiles=guests[guest],
                                       N=1, optimize=True)
        # save to mol file
        Chem.MolToMolFile(mol, filename=guest + '.mol')
        # save to xyz file
        Chem.MolToPDBFile(mol, filename=guest + '.pdb')
        convert_PDB_2_XYZ(file=guest + '.pdb')
