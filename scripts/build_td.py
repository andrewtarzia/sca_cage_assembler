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
            'amine_smiles': 'Nc1ccc(-c2ccc(N)cc2)cc1',
            'alde_smiles': 'O=Cc1ccccn1'
        },
        'lig2': {
            'amine_smiles': 'Nc1ccc(-c2ccc(-c3ccc(N)cc3)cc2)cc1',
            'alde_smiles': 'O=Cc1ccccn1'
        },
        'tripod': {
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


def build_complexes():
    metal = build_metal()
    n_atom = build_N_atom()
    metal_centre = build_metal_centre(metal, n_atom)
    metal_centre = stk.BuildingBlock.init_from_molecule(
        metal_centre,
        # functional_groups=['metal_bound_NMN']
        functional_groups=['metal_bound_N']
    )

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

    # print('doing rdkit optimisation')
    # optimizer = stk.MetalOptimizer(
    #     metal_binder_distance=2.0,
    #     metal_binder_fc=1.0e3,
    #     binder_ligand_fc=0,
    #     ignore_vdw=False,
    #     rel_distance=None,
    #     res_steps=100,
    #     restrict_bonds=True,
    #     restrict_angles=True,
    #     restrict_orientation=True,
    #     max_iterations=20,
    #     do_long_opt=False
    # )
    # optimizer.optimize(s_complex)
    # s_complex.write(f's_complex_rdkit.mol')
    # s_complex.write(f's_complex_rdkit.xyz')
    # s_complex.dump(f's_complex_rdkit.json')
    # sys.exit()
    print('doing UFF4MOF optimisation')
    gulp_opt = stk.GulpMetalOptimizer(
        gulp_path='/home/atarzia/software/gulp-5.1/Src/gulp/gulp',
        metal_FF='Fe6+2',
        output_dir=f's_complex_uff1'
    )
    gulp_opt.assign_FF(s_complex)
    gulp_opt.optimize(mol=s_complex)
    s_complex.write(f's_complex_uff4mof.mol')
    s_complex.write(f's_complex_uff4mof.xyz')
    s_complex.dump(f's_complex_uff4mof.json')

    S_complex = stk.BuildingBlock.init_from_molecule(
        s_complex,
        functional_groups=['bromine'],
    )
    S_complex.write('built_s.mol')
    print(S_complex)
    print(S_complex.func_groups)

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
    print(R_complex)
    print(R_complex.func_groups)
    R_complex.write('built_r.mol')

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
    return
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
    S_complex, R_complex = build_complexes()
    print(S_complex, R_complex)

    n_metals = [4]
    topologies = {
        'm4l6_td': stk.cage.M4L6_Oct(),
        # 'm4l4_td': (stk.cage.M4L4_Oct(), 4),
        # 'td4octssss': stk.cage.M4L6_Oct_SSSS(),
    }
    ratios =[(4, 0), (3, 1), (2, 2), (1, 3), (0, 4)]
    complexes = [S_complex, R_complex]
    tripod = ligands['tripod']['molecule']
    for rat in ratios:
        print(rat[0], rat[1], rat[0]+rat[1])
        for i, topo in enumerate(topologies):
            top = topologies[topo]
            cage_name = f"{rat[0]}_{rat[1]}_{topo}"
            if not exists(f'{cage_name}_opt.mol'):
                print(f'build {cage_name}')
                # Build homo leptic cages.
                cage = stk.ConstructedMolecule(
                    building_blocks=[
                        complexes[0],
                        complexes[1],
                        tripod
                    ],
                    topology_graph=top,
                    building_block_vertices={
                        complexes[0]: top.vertices[:rat[0]],
                        complexes[1]: top.vertices[
                            rat[0]:rat[0]+rat[1]
                        ],
                        tripod: top.vertices[n_metals[i]:]
                    }
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
            continue
            if not exists(f'{cage_name}_opt.ey'):
                calculate_energy(cage_name, n_metals=n_metals[i])


if __name__ == "__main__":
    main()
