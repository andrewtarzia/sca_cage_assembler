#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

import sys
from os.path import exists
import stk
from rdkit.Chem import AllChem as Chem


def build_ligand():
    amine = stk.BuildingBlock(
        smiles='Nc1ccc(-c2ccc(N)cc2)cc1',
        functional_groups=['amine']
    )
    aldehyde = stk.BuildingBlock(
        smiles='O=Cc1ccccn1',
        functional_groups=['aldehyde']
    )

    p_top = stk.polymer.Linear('ABA', 1)
    polymer = stk.ConstructedMolecule(
        building_blocks=[aldehyde, amine],
        topology_graph=p_top
    )
    opt = stk.UFF()
    opt.optimize(polymer)
    polymer.write('polymer.mol')
    lig_mol = stk.BuildingBlock.init_from_molecule(
        polymer,
        functional_groups=['CNC_metal']
        # functional_groups=['NCCN_metal']
    )

    ligand = {
        'molecule': lig_mol,
        'name': 'lig1',
    }

    return ligand


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


def build_homoleptic_cage(
    metal,
    ligand,
    cage_name,
    top,
    n_metals,
    metal_type
):
    cage = stk.ConstructedMolecule(
        building_blocks=[metal, ligand],
        topology_graph=top,
        building_block_vertices={
            metal: top.vertices[:n_metals],
            ligand: top.vertices[n_metals:],
        }
    )
    cage.write(f'{cage_name}_unopt.mol')
    cage.write(f'{cage_name}_unopt.xyz')
    cage.dump(f'{cage_name}_unopt.json')

    # optimizer = stk.MetalOptimizer(
    #     metal_binder_distance=mb_dist,
    #     metal_binder_fc=mb_fc,
    #     binder_ligand_fc=bl_fc,
    #     ignore_vdw=False,
    #     rel_distance=rd,
    #     res_steps=steps,
    #     restrict_bonds=True,
    #     restrict_angles=True,
    #     restrict_orientation=True,
    #     max_iterations=40,
    #     do_long_opt=do_long
    # )
    print('doing UFF4MOF optimisation')
    gulp_opt = stk.GulpMetalOptimizer(
        gulp_path='/home/atarzia/software/gulp-5.1/Src/gulp/gulp',
        metal_FF=metal_type,
        output_dir=f'cage_opt_{cage_name}_uff1'
    )
    gulp_opt.assign_FF(cage)
    gulp_opt.optimize(mol=cage)
    cage.write(f'{cage_name}_uff4mof.mol')
    cage.write(f'{cage_name}_uff4mof.xyz')
    cage.dump(f'{cage_name}_uff4mof.json')

    print('doing UFF4MOF optimisation 2')
    gulp_opt2 = stk.GulpMetalOptimizer(
        gulp_path='/home/atarzia/software/gulp-5.1/Src/gulp/gulp',
        metal_FF=metal_type,
        output_dir=f'cage_opt_{cage_name}_uff2'
    )
    gulp_opt2.assign_FF(cage)
    gulp_opt2.optimize(mol=cage)
    cage.write(f'{cage_name}_uff2.mol')
    cage.write(f'{cage_name}_uff2.xyz')
    cage.dump(f'{cage_name}_uff2.json')

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
        N_conformers=20
    )
    gulp_MD.assign_FF(cage)
    gulp_MD.optimize(cage)
    cage.write(f'{cage_name}_MD.mol')
    cage.write(f'{cage_name}_MD.xyz')
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
    cage.write(f'{cage_name}_prextb.xyz')
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
    cage.write(f'{cage_name}_opt.xyz')
    cage.dump(f'{cage_name}_opt.json')


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


def main():
    metal = build_metal()
    n_atom = build_N_atom()
    metal_centre = build_metal_centre(metal, n_atom)
    metal_centre = stk.BuildingBlock.init_from_molecule(
        metal_centre,
        # functional_groups=['metal_bound_NMN']
        functional_groups=['metal_bound_N']
    )

    metal_centre.write('metal_centre.mol')

    print(metal_centre.func_groups)

    ligand = build_ligand()
    print(ligand['molecule'].func_groups)

    n_metals = [4]
    topologies = {
        'td4oct': stk.cage.M4L6_Oct(),
    }
    for i, topo in enumerate(topologies):
        top = topologies[topo]
        cage_name = f"{ligand['name']}_{topo}"
        if not exists(f'{cage_name}_opt.mol'):
            print(f'build {cage_name}')
            # Build homo leptic cages.
            build_homoleptic_cage(
                metal=metal_centre,
                ligand=ligand['molecule'],
                cage_name=cage_name,
                top=top,
                n_metals=n_metals[i],
                metal_type='Fe6+2'
            )
        if not exists(f'{cage_name}_opt.ey'):
            calculate_energy(cage_name, n_metals=n_metals[i])


if __name__ == "__main__":
    main()
