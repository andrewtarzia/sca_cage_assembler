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


def available_topologies(string):
    """
    Get stk function of desired topology.

    """

    topologies = {
        'oct_lambda': stk.cage.Octahedral_Lambda(),
        'oct_delta': stk.cage.Octahedral_Delta()
    }

    try:
        return topologies[string]
    except KeyError:
        raise KeyError(f'{string} not in {topologies.keys()}')


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


def metal_FFs(CN):
    """
    Define metal FF names for UFF4MOF.

    Key = Atomic number
    Value = UFF4MOF type
    CN = coordination number of metal.

    """
    if CN == 4:
        dicts = {
            26: 'Fe4+2',
            27: 'Co4+2',
            28: 'Ni4+2',
            30: 'Zn4+2',
            42: 'Mo4f2',
            45: 'Rh6+3',  # No alternative available.
            46: 'Pd4+2',
            48: 'Cd4f2',
            78: 'Pt4+2',
        }
    elif CN == 6:
        dicts = {
            26: 'Fe6+2',
            27: 'Co6+2',
            28: 'Ni4+2',  # No alternative available.
            30: 'Zn4+2',  # No alternative available for 90 degrees.
            45: 'Rh6+3',
            48: 'Cd4f2',  # No alternative available.
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
        Stk BuildingBlock of metal atom.

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


def optimize_SCA_complex(complex, name, dict, metal_FFs):
    """
    Optimize a sub-component self assmebly complex.

    """

    print(f'doing UFF4MOF optimisation for {name}')
    gulp_opt = stk.GulpMetalOptimizer(
        gulp_path='/home/atarzia/software/gulp-5.1/Src/gulp/gulp',
        metal_FF=metal_FFs,
        output_dir=f'{name}_uff1'
    )
    gulp_opt.assign_FF(complex)
    gulp_opt.optimize(mol=complex)
    complex.write(f'{name}_uff1.mol')
    complex.dump(f'{name}_uff1.json')

    print(f'doing xTB optimisation for {name}')
    xtb_opt = stk.XTB(
        xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
        output_dir=f'{name}_xtb',
        gfn_version=2,
        num_cores=6,
        opt_level='tight',
        charge=dict['total_charge'],
        num_unpaired_electrons=dict['unpaired_e'],
        max_runs=1,
        calculate_hessian=False,
        unlimited_memory=True
    )
    xtb_opt.optimize(mol=complex)
    complex.write(f'{name}_opt.mol')
    complex.dump(f'{name}_opt.json')

    return complex


def build_SCA_complex(
    metal_centre,
    bidentate_ligand,
    complex_top,
    name,
):
    """
    Build an stk metal complex post SCA coordination.

    Parameters
    ----------
    metal_centre : :class:`stk.BuildingBlock`
        Stk BuildingBlock of metal centre.

    bidentate_ligand : :class:`stk.BuildingBlock`
        Ligand to bind to metal centre.

    complex_top : :class:`stk.Topology`
        Topology of metal complex to build.

    name : :class:`str`
        Name of structure.

    Returns
    -------
    complex : :class:`stk.BuildingBlock`
        Built stk molecule as :class:`stk.BuildingBlock` with
        `return_FG` FGs.

    """
    complex = stk.ConstructedMolecule(
        building_blocks=[
            metal_centre,
            bidentate_ligand
        ],
        topology_graph=complex_top,
        building_block_vertices={
            metal_centre: tuple([complex_top.vertices[0]]),
            bidentate_ligand: complex_top.vertices[1:]
        }
    )
    complex.write(f'{name}.mol')
    complex.write(f'{name}.mol')

    return complex
