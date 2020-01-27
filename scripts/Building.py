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


def experimentally_tested(mol_name):
    """
    A dataset of experimentally built cages.

    THIS NEEDS DEFINITION BASED ON CAGE NOMENCLATURE.

    Here a combo is defined as by all sub-component self-assembly
    products:
    i.e. one cage would be formed by a mixture of:
        (SCA, quadX) or (SCA, triX)

    """

    raise NotImplementedError()

    tested_ligand_combos = []

    return tested_ligand_combos


def metal_FFs():
    """
    Define metal FF names for UFF4MOF.

    Key = Atomic number
    Value = UFF4MOF type

    """
    dicts = {
        30: 'Zn4+2', 28: 'Ni4+2',
        78: 'Pt4+2', 46: 'Pd4+2',
        45: 'Rh6+3', 42: 'Mo4f2'
    }

    return dicts


def build_atom(smiles, FG):
    """
    Build an stk readable atom using RDKit.

    Parameters
    ----------
    smiles : :class:`str`
        Smiles of metal atom to use - should include charge.

    FG : :class:`str`
        Type of functional group to assign to atom.

    Returns
    -------
    stk_atom : :class:`stk.BuildingBlock`
        Built stk molecule with functional group.

    """

    atom = rdkit.MolFromSmiles(smiles)
    atom.AddConformer(rdkit.Conformer(atom.GetNumAtoms()))
    stk_atom = stk.BuildingBlock.init_from_rdkit_mol(
        atom,
        functional_groups=[FG],
    )
    return stk_atom


def build_metal(metal_smiles, no_fgs):
    """
    Build an stk readable metal atom using RDKit.

    Parameters
    ----------
    metal_smiles : :class:`str`
        Smiles of metal atom to use - should include charge.

    no_fgs : :class:`int`
        Number of functional groups to give metal atom.

    Returns
    -------
    metal : :class:`stk.BuildingBlock`
        Built stk molecule with functional groups.

    """

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


def build_metal_centre(metal, topology, binding_atom, return_FG):
    """
    Build an stk metal centre.

    Parameters
    ----------
    metal : :class:`stk.BuildingBlock`
        Smiles of metal atom to use - should include charge.

    topology : :class:`stk.MetalCentre`
        Topology of metal centre.

    binding_atom : :class:`stk.BuildingBlock`
        Atom to be placed at coordinating sites.

    return_FG : :class:`list`
        Functional groups to read in building block.


    Returns
    -------
    complex : :class:`stk.BuildingBlock`
        Built stk molecule as :class:`stk.BuildingBlock` with
        `return_FG` FGs.

    """

    complex = stk.ConstructedMolecule(
        building_blocks=[metal, binding_atom],
        topology_graph=topology,
        building_block_vertices={
            metal: tuple([topology.vertices[0]]),
            binding_atom: topology.vertices[1:]
        }
    )

    complex = stk.BuildingBlock.init_from_molecule(
        complex,
        functional_groups=return_FG
    )

    return complex
