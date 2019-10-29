#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

import sys
from os.path import exists
import stk
from rdkit.Chem import AllChem as Chem


def build_ligands():
    ligands = {
        'lig1': {
            'amine_smiles': (
                'Nc1ccc(-c2ccc(N)cc2S(=O)(=O)O)c(S(=O)(=O)O)c1'
            ),
            'alde_smiles': 'O=Cc1ccccn1'
        },
        'lig2': {
            'amine_smiles': 'Nc1ccc(-c2ccc(-c3ccc(N)cc3)cc2)cc1',
            'alde_smiles': 'O=Cc1ccccn1'
        },
        'spacer3': {'name': 's3'},
        'spacer1': {'name': 's1'},
        'spacerni': {'name': 'sNi'},
        'bident1': {'name': 'b1'},
        'bident2': {'name': 'b2'},
    }

    for ligand in ligands:
        if 'spacer' in ligand:
            if ligand == 'spacerni':
                continue
            lig_mol = stk.BuildingBlock.init_from_file(
                f'{ligand}.mol',
                functional_groups=['bromine']
            )
            ligands[ligand]['molecule'] = lig_mol
        elif 'bident' in ligand:
            lig_mol = stk.BuildingBlock.init_from_file(
                f'{ligand}.mol',
                functional_groups=['CNC_metal']
            )
            ligands[ligand]['molecule'] = lig_mol
        else:
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
            if 'lig' not in ligand:
                # spacer or tripod
                lig_mol = stk.BuildingBlock.init_from_molecule(
                    polymer,
                    functional_groups=['bromine']
                    # functional_groups=['NCCN_metal']
                )
            else:
                lig_mol = stk.BuildingBlock.init_from_molecule(
                    polymer,
                    functional_groups=['CNC_metal']
                    # functional_groups=['NCCN_metal']
                )
            ligands[ligand]['molecule'] = lig_mol
            lig_mol.write(ligand+'.mol')
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


def calculate_energy(cage_name, n_metals):
    if not exists(f'{cage_name}_opt.json'):
        sys.exit(f'{cage_name}_opt.json does not exist.')
    else:
        cage_file = f'{cage_name}_opt.json'
        cage = stk.ConstructedMolecule.load(cage_file)

    # Extract energy.
    xtb_energy = stk.XTBEnergy(
        xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
        output_dir=f'cage_opt_{cage_name}',
        num_cores=6,
        charge=n_metals*2,
        num_unpaired_electrons=0,
        electronic_temperature=300,
        unlimited_memory=True
    )
    energy = xtb_energy.get_energy(cage)

    # Save to .ey file.
    energy_file = f'{cage_name}_opt.ey'
    with open(energy_file, 'w') as f:
        f.write(f'{energy}\n')


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
    complex.write(f'{name}.mol')

    if exists(f'{name}_opt.mol'):
        complex.update_from_file(f'{name}_opt.mol')
    else:
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

        print('doing XTB optimisation')
        xtb_opt = stk.XTB(
            xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
            output_dir=f'{name}_xtb',
            gfn_version=2,
            num_cores=6,
            opt_level='tight',
            charge=2,
            num_unpaired_electrons=0,
            max_runs=1,
            electronic_temperature=1000,
            calculate_hessian=False,
            unlimited_memory=True
        )
        xtb_opt.optimize(mol=complex)
        complex.write(f'{name}_opt.mol')
        complex.write(f'{name}_opt.xyz')
        complex.dump(f'{name}_opt.json')

        # print('doing UFF4MOF optimisation 2')
        # gulp_opt3 = stk.GulpMetalOptimizer(
        #     gulp_path='/home/atarzia/software/gulp-5.1/Src/gulp/gulp',
        #     metal_FF='Fe6+2',
        #     output_dir=f's_complex_uff3'
        # )
        # gulp_opt3.assign_FF(complex)
        # gulp_opt3.optimize(mol=complex)
        # complex.write(f's_complex_postxtb.mol')
        # complex.write(f's_complex_4.xyz')
        # complex.dump(f's_complex_postxtb.json')

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
        # functional_groups=['metal_bound_NMN']
        functional_groups=['metal_bound_N']
    )

    metal_centre.write('metal_centre.mol')

    bidentate_ligand = bidentate['molecule']
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


def optimize_cage(cage, cage_name, n_metals, metal_type):

    # print('doing OPLS optimisation')
    # optimizer = stk.MacroModelFFMetalOptimizer(
    #     macromodel_path='/home/atarzia/software/schrodinger_install',
    #     output_dir=f'{cage_name}_opls',
    #     restrict_all_bonds=True
    # )
    # optimizer.optimize(cage)
    # cage.write(f'{cage_name}_opls.mol')
    # cage.write(f'{cage_name}_2.xyz')
    # cage.dump(f'{cage_name}_opls.json')
    #
    # print('doing OPLS MD optimisation')
    # optimizer = stk.MacroModelMDMetalOptimizer(
    #     macromodel_path='/home/atarzia/software/schrodinger_install',
    #     output_dir=f'{cage_name}_oplsMD',
    #     restrict_all_bonds=False
    # )
    # optimizer.optimize(cage)
    # cage.write(f'{cage_name}_oplsMD.mol')
    # cage.write(f'{cage_name}_3.xyz')
    # cage.dump(f'{cage_name}_oplsMD.json')
    #
    # import sys
    # sys.exit()

    # print('doing rdkit optimisation')
    # optimizer = stk.MetalOptimizer(
    #     metal_binder_distance=2.0,
    #     metal_binder_fc=1.0e2,
    #     binder_ligand_fc=0,
    #     ignore_vdw=False,
    #     rel_distance=None,
    #     res_steps=100,
    #     restrict_bonds=True,
    #     restrict_angles=True,
    #     restrict_orientation=True,
    #     max_iterations=40,
    #     do_long_opt=True
    # )
    # optimizer.optimize(cage)
    # cage.write(f'{cage_name}_rdkit.mol')
    # cage.write(f'{cage_name}_rdkit.xyz')
    # cage.dump(f'{cage_name}_rdkit.json')

    print('doing UFF4MOF optimisation')
    gulp_opt = stk.GulpMetalOptimizer(
        gulp_path='/home/atarzia/software/gulp-5.1/Src/gulp/gulp',
        metal_FF=metal_type,
        output_dir=f'cage_opt_{cage_name}_uff1'
    )
    gulp_opt.assign_FF(cage)
    gulp_opt.optimize(mol=cage)
    cage.write(f'{cage_name}_uff4mof.mol')
    cage.write(f'{cage_name}_2.xyz')
    cage.dump(f'{cage_name}_uff4mof.json')

    print('doing UFF4MOF MD')
    gulp_MD = stk.GulpMDMetalOptimizer(
        gulp_path='/home/atarzia/software/gulp-5.1/Src/gulp/gulp',
        metal_FF=metal_type,
        output_dir=f'cage_opt_{cage_name}_MD',
        integrator='stochastic',
        ensemble='nvt',
        temperature='700',
        equilbration='1.0',
        production='5.0',
        timestep='0.5',
        N_conformers=20,
        opt_conformers=False
    )
    gulp_MD.assign_FF(cage)
    gulp_MD.optimize(cage)
    cage.write(f'{cage_name}_MD.mol')
    cage.write(f'{cage_name}_3.xyz')
    cage.dump(f'{cage_name}_MD.json')

    print('doing UFF4MOF optimisation 3')
    gulp_opt3 = stk.GulpMetalOptimizer(
        gulp_path='/home/atarzia/software/gulp-5.1/Src/gulp/gulp',
        metal_FF=metal_type,
        output_dir=f'cage_opt_{cage_name}_uff3'
    )
    gulp_opt3.assign_FF(cage)
    gulp_opt3.optimize(mol=cage)
    cage.write(f'{cage_name}_prextb.mol')
    cage.write(f'{cage_name}_4.xyz')
    cage.dump(f'{cage_name}_prextb.json')

    print('doing XTB optimisation')
    xtb_opt = stk.XTB(
        xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
        output_dir=f'cage_opt_{cage_name}',
        gfn_version=2,
        num_cores=6,
        opt_level='tight',
        charge=n_metals*2,
        num_unpaired_electrons=0,
        max_runs=1,
        electronic_temperature=1000,
        calculate_hessian=False,
        unlimited_memory=True
    )
    xtb_opt.optimize(mol=cage)
    cage.write(f'{cage_name}_opt.mol')
    cage.write(f'{cage_name}_5.xyz')
    cage.dump(f'{cage_name}_opt.json')


def main():
    ligands = build_ligands()

    n_metals = [4, 4, 4]
    cages = {
        'm4l6s_td': (
            stk.cage.M4L6_Oct_Spacer(),
            ligands['bident2'],
            ligands['spacer1']
        ),
        'm4l6_td': (
            stk.cage.M4L6_Oct(),
            ligands['bident1'],
            None
        ),
        'm4l4_td': (
            stk.cage.M4L4_Oct_Spacer(),
            ligands['bident2'],
            ligands['spacer3']
        ),
    }
    ratios =[(4, 0), (3, 1), (2, 2), (1, 3), (0, 4)]
    for rat in ratios:
        print(f'ratio S:R: {rat[0]}:{rat[1]}')
        for i, topo in enumerate(cages):
            top, bident, spacer = cages[topo]
            bident_name = bident['name']
            S_complex, R_complex = build_complexes(bident)
            complexes = [S_complex, R_complex]
            if spacer is None:
                cage_name = f"{bident_name}_{rat[0]}_{rat[1]}_{topo}"
                bbs = [complexes[0], complexes[1]]
                bb_vs = {
                    complexes[0]: top.vertices[:rat[0]],
                    complexes[1]: top.vertices[rat[0]:rat[0]+rat[1]]
                }
            else:
                print(spacer)
                spacer_mol = spacer['molecule']
                s_name = spacer['name']
                cage_name = (
                    f"{bident_name}_{s_name}_{rat[0]}_{rat[1]}_{topo}"
                )
                bbs = [complexes[0], complexes[1], spacer_mol]
                bb_vs = {
                    complexes[0]: top.vertices[:rat[0]],
                    complexes[1]: top.vertices[rat[0]:rat[0]+rat[1]],
                    spacer_mol: top.vertices[n_metals[i]:]
                }
            if not exists(f'{cage_name}_opt.mol'):
                print(f'build {cage_name}')
                # Build homo leptic cages.
                cage = stk.ConstructedMolecule(
                    building_blocks=bbs,
                    topology_graph=top,
                    building_block_vertices=bb_vs
                )
                # print(cage.func_groups)
                cage.write(f'{cage_name}_unopt.mol')
                cage.write(f'{cage_name}_1.xyz')
                cage.dump(f'{cage_name}_unopt.json')
                optimize_cage(
                    cage,
                    cage_name=cage_name,
                    n_metals=n_metals[i],
                    metal_type='Fe6+2'
                )
            if not exists(f'{cage_name}_opt.ey'):
                calculate_energy(cage_name, n_metals=n_metals[i])


if __name__ == "__main__":
    main()
