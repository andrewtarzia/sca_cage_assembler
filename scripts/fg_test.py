#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

import stk
from rdkit.Chem import AllChem as Chem


def build_ligands():
    ligands = {
        'bident1': {'name': 'b1'},
    }

    for ligand in ligands:
        if 'bident' in ligand:
            lig_mol = stk.BuildingBlock.init_from_file(
                f'{ligand}.mol',
                functional_groups=['CNC_metal']
            )
            ligands[ligand]['molecule'] = lig_mol

    return ligands


def build_metal():
    m = Chem.MolFromSmiles('[Fe+2]')
    m.AddConformer(Chem.Conformer(m.GetNumAtoms()))
    metal = stk.BuildingBlock.init_from_rdkit_mol(
        m,
        functional_groups=None,
    )
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


def build_complex(metal_centre, bidentate_ligand, complex_top, name):
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
    complex.write(f'{name}.pdb')

    print(f'doing opt for {name}')
    print('doing UFF4MOF optimisation')
    gulp_opt = stk.GulpMetalOptimizer(
        gulp_path='/home/atarzia/software/gulp-5.1/Src/gulp/gulp',
        metal_FF='Fe6+2',
        output_dir=f'{name}_uff1'
    )
    gulp_opt.assign_FF(complex)
    gulp_opt.optimize(mol=complex)
    complex.write(f'{name}_uff1.mol')
    complex.write(f'{name}_uff1.xyz')
    complex.dump(f'{name}_uff1.json')

    # print('doing UFF4MOF optimisation')
    # gulp_opt = stk.GulpMetalOptimizer(
    #     gulp_path='/home/atarzia/software/gulp-5.1/Src/gulp/gulp',
    #     metal_FF='Fe6+2',
    #     output_dir=f'{name}_uff1'
    # )
    # gulp_opt.assign_FF(complex)
    # gulp_opt.optimize(mol=complex)
    # complex.write(f'{name}_uff2.mol')
    # complex.write(f'{name}_uff2.xyz')
    # complex.dump(f'{name}_uff2.json')

    # print('doing XTB optimisation')
    # xtb_opt = stk.XTB(
    #     xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
    #     output_dir=f'{name}_xtb',
    #     gfn_version=2,
    #     num_cores=6,
    #     opt_level='tight',
    #     charge=2,
    #     num_unpaired_electrons=0,
    #     max_runs=1,
    #     electronic_temperature=1000,
    #     calculate_hessian=False,
    #     unlimited_memory=True
    # )
    # xtb_opt.optimize(mol=complex)
    # complex.write(f'{name}_xtb.mol')
    # complex.write(f'{name}_xtb.xyz')
    # complex.dump(f'{name}_xtb.json')

    complex = stk.BuildingBlock.init_from_molecule(
        complex,
        functional_groups=['bromine'],
    )

    return complex


def build_complexes(bidentate):
    print('building complexes')
    metal = build_metal()
    n_atom = build_N_atom()
    metal_centre = build_metal_centre(metal, n_atom)
    metal_centre = stk.BuildingBlock.init_from_molecule(
        metal_centre,
        functional_groups=['metal_bound_N']
    )
    # new_fgs = []
    # bonders = []
    # for i in metal_centre.func_groups:
    #     if i.bonders[0] not in bonders and i.bonders[1] not in bonders:
    #         new_fgs.append(i)
    #         bonders.append(i.bonders[0])
    #         bonders.append(i.bonders[1])
    # metal_centre.func_groups = tuple(new_fgs)
    # print(metal_centre.func_groups)

    metal_centre.write('metal_centre.mol')

    bidentate_ligand = bidentate['molecule']
    print('------------')
    print(bidentate_ligand.func_groups)
    print('------------')
    print('------------')
    print(metal_centre.func_groups)
    print('------------')
    bident_name = bidentate['name']
    topologies = {
        's_c': stk.cage.Octahedral_S(),
        'r_c': stk.cage.Octahedral_R(),
    }

    for top in topologies:
        complex_top = topologies[top]
        comp = build_complex(
            metal_centre,
            bidentate_ligand,
            complex_top,
            name=f'{top}_{bident_name}'
        )
        if top == 's_c':
            S_complex = comp
        elif top == 'r_c':
            R_complex = comp

    print('done')
    return S_complex, R_complex


def main():
    ligands = build_ligands()
    bident = ligands['bident1']
    S_complex, R_complex = build_complexes(bident)


if __name__ == "__main__":
    main()
