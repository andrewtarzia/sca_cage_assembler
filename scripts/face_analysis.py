#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build and analyse faces from ligand library.

Author: Andrew Tarzia

Date Created: 19 Oct 2020

"""

import sys
import json
from os.path import exists, join
from itertools import combinations
import matplotlib.pyplot as plt
from glob import glob
import numpy as np

import stk
import stko

from molecule_building import metal_FFs
from cubeface import CubeFace

from atools import (
    MOC_collapse_mc,
    MOC_uff_opt,
    get_atom_distance,
    colors_i_like,
)


def load_complex(filename):

    fgfactory = stk.SmartsFunctionalGroupFactory(
        smarts='[#7X3]~[#6]~[#6]~[#7X3]~[#35]',
        bonders=(3, ),
        deleters=(4, ),
        placers=(0, 1, 2, 3),
    )

    name = filename.replace('.mol', '')
    # Need to define more than one placer id for complexes to ensure
    # alignment -- use the NCCN plane (attached to the Br) to define
    # the orientation of the complex.
    complex = stk.BuildingBlock.init_from_file(
        filename,
        functional_groups=[fgfactory]
    )

    return name, complex


def optimize_complex(complex, name):

    opt_name = f'{name}_opt.mol'
    if exists(opt_name):
        return complex.with_structure_from_file(opt_name)
    else:
        print(f'doing UFF4MOF optimisation for {name}')
        gulp_opt = stko.GulpUFFOptimizer(
            gulp_path='/home/atarzia/software/gulp-5.1/Src/gulp/gulp',
            metal_FF=metal_FFs(CN=6),
            output_dir=f'{name}_uff1'
        )
        gulp_opt.assign_FF(complex)
        complex = gulp_opt.optimize(mol=complex)
        complex.write(f'{name}_uff1.mol')

        print(f'doing xTB optimisation for {name}')
        xtb_opt = stko.XTB(
            xtb_path='/home/atarzia/software/xtb-6.3.1/bin/xtb',
            output_dir=f'{name}_xtb',
            gfn_version=2,
            num_cores=6,
            opt_level='tight',
            charge=2,
            num_unpaired_electrons=0,
            max_runs=1,
            calculate_hessian=False,
            unlimited_memory=True
        )
        complex = xtb_opt.optimize(mol=complex)
        complex.write(f'{name}_opt.mol')


def load_ligands(directory):

    ligands = {}
    for lig in glob(join(directory, 'quad2*_opt.mol')):
        l_name = lig.replace(directory, '').replace('_opt.mol', '')
        ligands[l_name] = stk.BuildingBlock.init_from_file(
            lig,
            functional_groups=[stk.BromoFactory()],
        )

    return ligands


def build_face(
    face_name,
    lig_structure,
    del_complex,
    lam_complex,
    face_topo
):

    face_file = f'{face_name}.mol'

    face = stk.ConstructedMolecule(
        topology_graph=CubeFace(
            building_blocks={
                del_complex: face_topo['d_pos'],
                lam_complex: face_topo['l_pos'],
                lig_structure: (4, ),
            },
            vertex_alignments=face_topo['va'],
        )
    )

    face.write(face_file)
    return face


def optimize_face(face, face_name):

    coll_file = f'{face_name}_coll.mol'
    opt_file = f'{face_name}_opt.mol'

    if exists(opt_file):
        return face.with_structure_from_file(opt_file)

    # Collapser MC algorithm.
    if exists(coll_file):
        opt_face = face.with_structure_from_file(coll_file)
    else:
        target_bond_length = 1.2
        num_steps = 2000
        step_size = 0.25
        opt_face = MOC_collapse_mc(
            cage=face,
            cage_name=face_name,
            step_size=step_size,
            target_bond_length=target_bond_length,
            num_steps=num_steps,
        )
        opt_face.write(coll_file)

    # Short restrained UFF opt.
    custom_metal_FFs = metal_FFs(CN=6)
    opt_face = MOC_uff_opt(
        opt_face,
        face_name,
        metal_FFs=custom_metal_FFs,
        CG=True,
        maxcyc=50,
        metal_ligand_bond_order='',
    )
    opt_face.write(opt_file)

    return opt_face


def get_all_bond_lengths(face):

    all_bls = []
    for bond in face.get_bonds():
        a1 = bond.get_atom1().get_id()
        a2 = bond.get_atom2().get_id()
        all_bls.append(get_atom_distance(
            molecule=face,
            atom1_id=a1,
            atom2_id=a2,
        ))

    return all_bls

    """
    Calculate geoemtrical properties of a face.

    """

    # Get metal atom ids.
    metal_atom_ids = []
    for atom in face.get_atoms():
        if atom.get_atomic_number() == metal_atomic_number:
            metal_atom_ids.append(atom.get_id())

    # Get metal-metal distances - as neighbours.
    metal_metal_distances = sorted(
        [(
            i,
            j,
            get_atom_distance(
                molecule=face,
                atom1_id=i,
                atom2_id=j
            ),
        ) for i, j in combinations(metal_atom_ids, r=2)
        ],
        key=lambda a: a[2]
    )
    properties = {'mismatches': []}
    for metal_atom_id in metal_atom_ids:
        # Assume neighbour is two shortest distances.
        # Save to dictionary -- atom_id: [dist1, dist2]
        neigh_dists = [
            i[2]
            for i in metal_metal_distances if metal_atom_id in i[:2]
        ][:2]
        properties[metal_atom_id] = neigh_dists
        mismatch = (
            100 * ((neigh_dists[1] - neigh_dists[0]) / neigh_dists[0])
        )
        properties['mismatches'].append(mismatch)

    return properties


def get_face_properties(face, face_name):

    json_file = f'{face_name}_properties.json'
    if exists(json_file):
        with open(json_file, 'r') as f:
            data = json.load(f)
    else:
        data = calculate_face_properties(face)
        with open(json_file, 'w') as f:
            json.dump(data, f)

    return data


def plot_face_mismatches(data, name):

    avg_mismatches = {}

    for face in data:
        avg_mismatch = np.average(data[face]['mismatches'])
        avg_mismatches[face] = avg_mismatch

    x_ticks = [int(i) for i in data]
    x_ticklabels = [i for i in data]

    width = 0.9
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.bar(
        x=[int(i) for i in avg_mismatches],
        height=[avg_mismatches[i] for i in avg_mismatches],
        width=width,
        facecolor=colors_i_like()[4],
        edgecolor='none',
        alpha=1,
    )

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylabel('face-side mismatch [%]', fontsize=16)
    ax.set_ylim(0, 200)
    # Set number of ticks for x-axis
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_ticklabels)

    fig.tight_layout()
    fig.savefig(
        f'f_mismatch_{name}.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def main():
    first_line = (
        'Usage: face_analysis.py '
        'lig_directory complex_name'
    )
    if (not len(sys.argv) == 3):
        print(f"""
{first_line}

    ligand_directory : (str)
        Directory with required ligand structures.

    complex_name : (str)
        Prefix of complex used in this run - defines face file names.

    """)
        sys.exit()
    else:
        ligand_directory = sys.argv[1]
        complex_name = sys.argv[2]

    # Define two metal building blocks (lambda, delta).
    complex_name = 'cl1'
    del_name, del_complex = load_complex(
        f'{complex_name}_zn_oct_del_face.mol'
    )
    lam_name, lam_complex = load_complex(
        f'{complex_name}_zn_oct_lam_face.mol'
    )

    # Optimise both complexes.
    del_complex = optimize_complex(del_complex, del_name)
    lam_complex = optimize_complex(lam_complex, lam_name)

    # Load in each ligand structure.
    ligands = load_ligands(ligand_directory)

    face_topologies = {
        '1': {
            'va': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
            'ratio': (2, 2),
            'd_pos': (1, 3),
            'l_pos': (0, 2),
        },
        '2': {
            'va': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
            'ratio': (2, 2),
            'd_pos': (0, 2),
            'l_pos': (1, 3),
        },
        '3': {
            'va': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
            'ratio': (4, 0),
            'd_pos': (),
            'l_pos': (0, 1, 2, 3),
        },
        '4': {
            'va': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
            'ratio': (3, 1),
            'd_pos': (3, ),
            'l_pos': (0, 1, 2),
        },
        '5': {
            'va': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
            'ratio': (3, 1),
            'd_pos': (2, ),
            'l_pos': (0, 1, 3, ),
        },
    }

    # Skip ligands not in database of 10.
    dataset_of_10 = [
        'quad2_1', 'quad2_12', 'quad2_2', 'quad2_3', 'quad2_8',
        'quad2_9', 'quad2_10', 'quad2_5', 'quad2_16', 'quad2_17',
    ]

    # Build and optimise five face options per ligand.
    for lig in sorted(ligands):
        if lig not in dataset_of_10:
            continue
        print(f'doing {lig}...')
        lig_structure = ligands[lig]
        # Build each face topology.
        lig_faces = {}
        for face_t in face_topologies:
            final_topology_dict = face_topologies[face_t]
            face_name = f'F_{complex_name}_{lig}_{face_t}'
            face = build_face(
                face_name=face_name,
                lig_structure=lig_structure,
                del_complex=del_complex,
                lam_complex=lam_complex,
                face_topo=final_topology_dict,
            )
            opt_face = optimize_face(face, face_name)
            all_bls = get_all_bond_lengths(opt_face)
            if max(all_bls) > 3:
                raise ValueError(
                    f'max bond length of {face_name} is {max(all_bls)}'
                )


            # Measure properties.
            face_properties = get_face_properties(opt_face, face_name)
            print(
                f":: {face_name}: "
                f"{np.average(face_properties['mismatches'])}"
            )
            lig_faces[face_t] = face_properties
        plot_face_mismatches(data=lig_faces, name=lig)

    sys.exit()


if __name__ == '__main__':
    main()
