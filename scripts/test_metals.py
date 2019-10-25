#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

import stk
from rdkit.Chem import AllChem as Chem


def build_ligands():
    ligands = {
        'lig1': {
            'amine_smiles': 'Nc1ccc(-c2ccc(N)cc2)cc1',
            'alde_smiles': 'O=Cc1ccccn1'
        },
        'lig2': {
            'amine_smiles': 'Nc1ccc(-c2ccc(-c3ccc(N)cc3)cc2)cc1',
            'alde_smiles': 'O=Cc1ccccn1'
        }
    }

    for ligand in ligands:
        li = ligands[ligand]
        amine = stk.BuildingBlock(
            smiles=li['amine_smiles'],
            # functional_groups=['amine_metal']
            functional_groups=['amine']
        )
        aldehyde = stk.BuildingBlock(
            smiles=li['alde_smiles'],
            # functional_groups=['pyridine_N_metal', 'aldehyde']
            functional_groups=['aldehyde']
        )
        amine.write('ami.mol')
        aldehyde.write('alde.mol')

        p_top = stk.polymer.Linear('ABA', 1)
        polymer = stk.ConstructedMolecule(
            building_blocks=[aldehyde, amine],
            topology_graph=p_top
        )
        opt = stk.UFF()
        opt.optimize(polymer)
        polymer.write(f'{ligand}.mol')
        lig_mol = stk.BuildingBlock.init_from_molecule(
            polymer,
            functional_groups=['CNC_metal']
            # functional_groups=['NCCN_metal']
        )
        ligands[ligand]['molecule'] = lig_mol
        ligands[ligand]['amine'] = amine
        ligands[ligand]['aldehyde'] = aldehyde

    return ligands


def build_metal():
    m = Chem.MolFromSmiles('[Fe+2]')
    m.AddConformer(Chem.Conformer(m.GetNumAtoms()))
    metal = stk.BuildingBlock.init_from_rdkit_mol(
        m,
        functional_groups=None,
    )
    print('---')
    for bond in metal.bonds:
        print('m bonds:', bond)
    print('---')
    metal_coord_info = {
        0: {
            'atom_ids': [0],
            'bonder_ids': [0],
            'deleter_ids': [None]
        },
        1: {
            'atom_ids': [0],
            'bonder_ids': [0],
            'deleter_ids': [None]
        },
        2: {
            'atom_ids': [0],
            'bonder_ids': [0],
            'deleter_ids': [None]
        },
        3: {
            'atom_ids': [0],
            'bonder_ids': [0],
            'deleter_ids': [None]
        },
        4: {
            'atom_ids': [0],
            'bonder_ids': [0],
            'deleter_ids': [None]
        },
        5: {
            'atom_ids': [0],
            'bonder_ids': [0],
            'deleter_ids': [None]
        },
    }
    metal = stk.assign_metal_fgs(
        building_block=metal,
        coordination_info=metal_coord_info
    )
    print('---')
    print('m', metal)
    print('---')
    return metal


def build_N_atom():
    m = Chem.MolFromSmiles('N')
    m.AddConformer(Chem.Conformer(m.GetNumAtoms()))
    n_atom = stk.BuildingBlock.init_from_rdkit_mol(
        m,
        functional_groups=['metal_bound_N'],
    )
    return n_atom


def build_metal_centre(metal, n_atom):
    m_top = stk.metal_centre.Octahedral()
    complex = stk.ConstructedMolecule(
        building_blocks=[metal, n_atom],
        topology_graph=m_top,
        building_block_vertices={
            metal: tuple([m_top.vertices[0]]),
            n_atom: m_top.vertices[1:]
        }
    )
    return complex


def build_complexes():
    metal = build_metal()
    n_atom = build_N_atom()
    metal_centre = build_metal_centre(metal, n_atom)
    print('---')
    print('mc', metal_centre)
    for bond in metal_centre.bonds:
        print('mcb', bond)
    metal_centre = stk.BuildingBlock.init_from_molecule(
        metal_centre,
        # functional_groups=['metal_bound_NMN']
        functional_groups=['metal_bound_N']
    )
    print('mc', metal_centre)
    for bond in metal_centre.bonds:
        print('mcb', bond)

    print('---')
    metal_centre.write('metal_centre.mol')
    bidentate_ligand = stk.BuildingBlock(
        'Brc1ccc(/N=C/c2ccccn2)cc1',
        functional_groups=['CNC_metal']
    )

    complex_top = stk.cage.Octahedral_S()
    s_complex = stk.ConstructedMolecule(
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
    s_complex.write('s_complex.mol')
    s_complex.write('s_complex.pdb')
    print('---')
    print('s', s_complex)
    for bond in s_complex.bonds:
        print('s', bond)
    S_complex = stk.BuildingBlock.init_from_molecule(
        s_complex,
        functional_groups=['bromine'],
    )
    S_complex.write('built_s.mol')
    print('S', S_complex)
    for bond in S_complex.bonds:
        print('S', bond)

    print('---')
    import sys
    sys.exit()
    # fg_makers = (stk.fg_types[name] for name in ['bromine'])
    # print(s_complex.func_groups)
    # s_complex_no_metals = s_complex.to_rdkit_mol_no_metals()
    # s_complex.func_groups = tuple(
    #     func_group
    #     for fg_maker in fg_makers
    #     for func_group in fg_maker.get_functional_groups(
    #         s_complex_no_metals
    #     )
    # )
    # print(s_complex.func_groups)

    complex_top = stk.cage.Octahedral_R()
    r_complex = stk.ConstructedMolecule(
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
    r_complex.write('r_complex.mol')
    r_complex.write('r_complex.pdb')

    print('doing UFF4MOF optimisation')
    gulp_opt = stk.GulpMetalOptimizer(
        gulp_path='/home/atarzia/software/gulp-5.1/Src/gulp/gulp',
        metal_FF='Fe6+2',
        output_dir=f'r_complex_uff1'
    )
    gulp_opt.assign_FF(r_complex)
    gulp_opt.optimize(mol=r_complex)
    r_complex.write(f'r_complex_uff4mof.mol')
    r_complex.write(f'r_complex_uff4mof.xyz')
    r_complex.dump(f'r_complex_uff4mof.json')

    R_complex = stk.BuildingBlock.init_from_molecule(
        r_complex,
        functional_groups=['bromine'],
    )

    return S_complex, R_complex


def main():
    ligands = build_ligands()
    S_complex, R_complex = build_complexes()


if __name__ == "__main__":
    main()
