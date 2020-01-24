#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Modules/functions for building structures.

Author: Andrew Tarzia

Date Created: 23 Jan 2020
"""

from rdkit.Chem import AllChem as rdkit
import stk


def build_metal(metal_smiles, no_fgs):
    m = rdkit.MolFromSmiles(metal_smiles)
    m.AddConformer(rdkit.Conformer(m.GetNumAtoms()))
    metal = stk.BuildingBlock.init_from_rdkit_mol(
        m,
        functional_groups=None,
    )
    fg_dict = {
        'atom_ids': [0],
        'bonder_ids': [0],
        'deleter_ids': [None]
    }
    metal_coord_info = {i: fg_dict for i in range(no_fgs)}

    metal = stk.assign_metal_fgs(
        building_block=metal,
        coordination_info=metal_coord_info
    )
    return metal
