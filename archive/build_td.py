#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

import sys
import itertools
from os.path import exists
import stk
from rdkit.Chem import AllChem as Chem
import env_set


def build_ligands():
    ligands = {
        'spacer3': {
            'name': 's3',
            'molecule': stk.BuildingBlock.init_from_file(
                'spacer3.mol',
                functional_groups=['bromine']
            )
        },
        'spacer1': {
            'name': 's1',
            'molecule': stk.BuildingBlock.init_from_file(
                'spacer1.mol',
                functional_groups=['bromine']
            )
        },
        'spacer1ni': {
            'name': '1ni',
            'molecule': stk.BuildingBlock.init_from_molecule(
                build_porphyrin(
                    metal='[Ni+2]',
                    porph_name='spacer1',
                    file_name='spacer1ni'
                ),
                functional_groups=['bromine']
            )
        },
        'bident1': {
            'name': 'b1',
            'molecule': stk.BuildingBlock.init_from_file(
                'bident1.mol',
                functional_groups=['CNC_metal']
            )
        },
        'bident2': {
            'name': 'b2',
            'molecule': stk.BuildingBlock.init_from_file(
                'bident2.mol',
                functional_groups=['CNC_metal']
            )
        }
    }

    return ligands


def build_porphyrin(metal, porph_name, file_name):
    metal = build_metal(metal_smiles=metal, no_fgs=4)
    porphyrin = stk.BuildingBlock.init_from_file(
        f'{porph_name}.mol',
        functional_groups=['metal_bound_N']
    )

    m_top = stk.metal_centre.Porphyrin()
    porph_m = stk.ConstructedMolecule(
        building_blocks=[metal, porphyrin],
        topology_graph=m_top,
        building_block_vertices={
            metal: tuple([m_top.vertices[0]]),
            porphyrin: m_top.vertices[1:]
        }
    )
    porph_m.write(f'{file_name}.mol')
    if exists(f'{file_name}_opt.mol'):
        porph_m.update_from_file(f'{file_name}_opt.mol')
    else:
        print('doing XTB optimisation')
        xtb_opt = stk.XTB(
            xtb_path=env_set.xtb_path(),
            output_dir=f'{file_name}',
            gfn_version=2,
            num_cores=6,
            opt_level='tight',
            charge=0,
            num_unpaired_electrons=0,
            max_runs=1,
            electronic_temperature=1000,
            calculate_hessian=False,
            unlimited_memory=True
        )
        xtb_opt.optimize(mol=porph_m)
        porph_m.write(f'{file_name}_opt.mol')

    return porph_m


def build_metal(metal_smiles, no_fgs):
    m = Chem.MolFromSmiles(metal_smiles)
    m.AddConformer(Chem.Conformer(m.GetNumAtoms()))
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
        xtb_path=env_set.xtb_path(),
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

    if exists(f'{name}_opt.mol'):
        complex.update_from_file(f'{name}_opt.mol')
    else:
        print(f'doing opt for {name}')
        print('doing UFF4MOF optimisation')
        gulp_opt = stk.GulpUFFOptimizer(
            gulp_path=env_set.gulp_path(),
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
            xtb_path=env_set.xtb_path(),
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
        # gulp_opt3 = stk.GulpUFFOptimizer(
        #     gulp_path=env_set.gulp_path(),
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

    metal = build_metal(metal_smiles='[Fe+2]', no_fgs=6)
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
        's_c': stk.cage.Octahedral_Lambda(),
        'r_c': stk.cage.Octahedral_Delta(),
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

    return S_complex, R_complex


def optimize_cage(cage, cage_name, n_metals, metal_types):

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
    gulp_opt = stk.GulpUFFOptimizer(
        gulp_path=env_set.gulp_path(),
        metal_FF=metal_types,
        output_dir=f'cage_opt_{cage_name}_uff1'
    )
    gulp_opt.assign_FF(cage)
    gulp_opt.optimize(mol=cage)
    cage.write(f'{cage_name}_uff4mof.mol')
    cage.write(f'{cage_name}_2.xyz')
    cage.dump(f'{cage_name}_uff4mof.json')

    print('doing UFF4MOF MD')
    gulp_MD = stk.GulpUFFMDOptimizer(
        gulp_path=env_set.gulp_path(),
        metal_FF=metal_types,
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
    gulp_opt3 = stk.GulpUFFOptimizer(
        gulp_path=env_set.gulp_path(),
        metal_FF=metal_types,
        output_dir=f'cage_opt_{cage_name}_uff3'
    )
    gulp_opt3.assign_FF(cage)
    gulp_opt3.optimize(mol=cage)
    cage.write(f'{cage_name}_prextb.mol')
    cage.write(f'{cage_name}_4.xyz')
    cage.dump(f'{cage_name}_prextb.json')

    print('doing XTB optimisation')
    xtb_opt = stk.XTB(
        xtb_path=env_set.xtb_path(),
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


def get_ratios(n_metals):
    rng = range(0, n_metals+1)
    rats = []
    for i in itertools.product(rng, rng):
        if i[0]+i[1] == n_metals:
            rats.append(i)
    return rats


def main():
    ligands = build_ligands()

    n_metals = [4, 4, 4, 4]
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
        'm4l6sni_td': (
            stk.cage.M4L6_Oct_Spacer(),
            ligands['bident2'],
            ligands['spacer1ni']
        ),
    }
    ratios = get_ratios(n_metals=4)
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
                    metal_types={26: 'Fe6+2', 28: 'Ni4+2'}
                )
            if not exists(f'{cage_name}_opt.ey'):
                calculate_energy(cage_name, n_metals=n_metals[i])


if __name__ == "__main__":
    main()