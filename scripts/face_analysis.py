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
import matplotlib.pyplot as plt
from glob import glob
import numpy as np

import stk
import stko

from molecule_building import metal_FFs
from cubeface import CubeFace
from facebuildingblock import FaceBuildingBlock, face_topology_dict
from utilities import get_query_atom_ids, get_atom_distance
import env_set


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
    metal_complex = stk.BuildingBlock.init_from_file(
        filename,
        functional_groups=[fgfactory]
    )

    return name, metal_complex


def optimize_complex(metal_complex, name):

    opt_name = f'{name}_opt.mol'
    if exists(opt_name):
        return metal_complex.with_structure_from_file(opt_name)
    else:
        print(f'doing UFF4MOF optimisation for {name}')
        gulp_opt = stko.GulpUFFOptimizer(
            gulp_path=env_set.gulp_path(),
            metal_FF=metal_FFs(CN=6),
            output_dir=f'{name}_uff1'
        )
        gulp_opt.assign_FF(metal_complex)
        metal_complex = gulp_opt.optimize(mol=metal_complex)
        metal_complex.write(f'{name}_uff1.mol')

        print(f'doing xTB optimisation for {name}')
        xtb_opt = stko.XTB(
            xtb_path=env_set.xtb_path(),
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
        metal_complex = xtb_opt.optimize(mol=metal_complex)
        metal_complex.write(f'{name}_opt.mol')
        return metal_complex


def load_ligands(directory):

    ligands = {}
    for lig in glob(join(directory, '*_opt.mol')):
        l_name = lig.replace(directory, '').replace('_opt.mol', '')
        bb = FaceBuildingBlock.init_from_file(
            lig,
            functional_groups=[stk.BromoFactory()],
        )
        if bb.get_num_functional_groups() == 4:
            ligands[l_name] = bb

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

        print(f'..doing collapser optimisation of {face_name}')
        output_dir = f'cage_opt_{face_name}_coll'
        optimizer = stko.CollapserMC(
            output_dir=output_dir,
            step_size=step_size,
            target_bond_length=target_bond_length,
            num_steps=num_steps,
        )
        opt_face = optimizer.optimize(mol=face)
        opt_face = opt_face.with_centroid([0, 0, 0])
        opt_face.write(coll_file)

    # Short restrained UFF opt.
    custom_metal_FFs = metal_FFs(CN=6)
    output_dir = f'cage_opt_{face_name}_uffCG'
    print(f'..doing UFF4MOF optimisation of {face_name}')
    gulp_opt = stko.GulpUFFOptimizer(
        gulp_path=env_set.gulp_path(),
        maxcyc=50,
        metal_FF=custom_metal_FFs,
        metal_ligand_bond_order='',
        output_dir=output_dir,
        conjugate_gradient=True
    )
    gulp_opt.assign_FF(opt_face)
    opt_face = gulp_opt.optimize(mol=opt_face)
    opt_face = opt_face.with_centroid([0, 0, 0])
    opt_face.write(opt_file)

    return opt_face


def long_optimize_face(face, face_name):
    gulp_file = f'{face_name}_lgulp.mol'
    lopt_file = f'{face_name}_lopt.mol'

    if exists(lopt_file):
        return face.with_structure_from_file(lopt_file)

    # Unrestrained UFF opt.
    custom_metal_FFs = metal_FFs(CN=6)
    output_dir = f'cage_opt_{face_name}_gulp'
    print(f'..doing UFF4MOF optimisation of {face_name}')
    gulp_opt = stko.GulpUFFOptimizer(
        gulp_path=env_set.gulp_path(),
        maxcyc=2000,
        metal_FF=custom_metal_FFs,
        metal_ligand_bond_order='',
        output_dir=output_dir,
        conjugate_gradient=False,
    )
    gulp_opt.assign_FF(face)
    opt_face = gulp_opt.optimize(mol=face)
    opt_face = opt_face.with_centroid([0, 0, 0])
    opt_face.write(gulp_file)

    # xTB opt.
    print(f'..doing XTB optimisation of {face_name}')
    xtb_opt = stko.XTB(
        xtb_path=env_set.xtb_path(),
        output_dir=f'cage_opt_{face_name}_xtb',
        gfn_version=2,
        num_cores=6,
        opt_level='crude',
        charge=8,
        num_unpaired_electrons=0,
        max_runs=1,
        electronic_temperature=300,
        calculate_hessian=False,
        unlimited_memory=True,
        solvent=None,
    )
    opt_face = xtb_opt.optimize(mol=opt_face)
    opt_face.write(lopt_file)

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


def calculate_face_properties(face, paths, face_type):
    """
    Calculate geometrical properties of a face.

    """

    properties = {
        'metals': None,
        'Ns': None,
        'Cs': None,
    }

    # Calculate N-N distance for this face based on long axis.
    path_atom_ids = paths['Ns']
    path1 = (
        (path_atom_ids[0], path_atom_ids[1]),
        (path_atom_ids[0], path_atom_ids[3]),
    )
    path2 = (
        (path_atom_ids[2], path_atom_ids[1]),
        (path_atom_ids[2], path_atom_ids[3]),
    )

    p1a_d = get_atom_distance(face, path1[0][0], path1[0][1])
    p1b_d = get_atom_distance(face, path1[1][0], path1[1][1])
    p2a_d = get_atom_distance(face, path2[0][0], path2[0][1])
    p2b_d = get_atom_distance(face, path2[1][0], path2[1][1])
    difference1 = abs(p1a_d-p1b_d)
    difference2 = abs(p2a_d-p2b_d)
    aspect_differences = (difference1, difference2)

    for path in paths:
        path_atom_ids = paths[path]
        mismatch_values = calculate_path_mismatch(
            face=face,
            path_ids=path_atom_ids,
            face_type=face_type,
        )

        properties[path] = {
            'aspect_differences': aspect_differences,
            'mms': (
                mismatch_values['mismatch1'],
                mismatch_values['mismatch2'],
            ),
            'dif': (
                mismatch_values['difference1'],
                mismatch_values['difference2'],
            ),
            'distances': {
                'path1': (
                    mismatch_values['p1a_d'],
                    mismatch_values['p1b_d'],
                ),
                'path2': (
                    mismatch_values['p2a_d'],
                    mismatch_values['p2b_d'],
                ),
            },
        }

    return properties


def show_long_axis(face, face_name):

    out_file = f'{face_name}_lashown.xyz'
    string = stk.XyzWriter().to_string(face)

    x_pos = np.linspace(-10, 10, 100)
    y_pos = [0 for i in x_pos]

    for x, y in zip(x_pos, y_pos):
        string += f'Ar {x} {y} 0\n'

    string = string.split('\n')
    string[0] = str(int(string[0]) + len(x_pos))
    string = '\n'.join(string)

    with open(out_file, 'w') as f:
        f.write(string)


def calculate_path_mismatch(face, path_ids, face_type):
    if face_type == 'i':
        path1a = (path_ids[2], path_ids[3])
        path1b = (path_ids[3], path_ids[0])
        path2a = (path_ids[0], path_ids[1])
        path2b = (path_ids[1], path_ids[2])
    elif face_type == 'ii':
        path1a = (path_ids[2], path_ids[3])
        path1b = (path_ids[3], path_ids[0])
        path2a = (path_ids[0], path_ids[1])
        path2b = (path_ids[1], path_ids[2])
    elif face_type == 'iii':
        path1a = (path_ids[3], path_ids[0])
        path1b = (path_ids[1], path_ids[0])
        path2a = (path_ids[2], path_ids[3])
        path2b = (path_ids[1], path_ids[2])
    elif face_type == 'iv':
        path1a = (path_ids[2], path_ids[3])
        path1b = (path_ids[3], path_ids[0])
        path2a = (path_ids[0], path_ids[1])
        path2b = (path_ids[1], path_ids[2])
    elif face_type == 'v':
        path1a = (path_ids[0], path_ids[1])
        path1b = (path_ids[1], path_ids[2])
        path2a = (path_ids[2], path_ids[3])
        path2b = (path_ids[3], path_ids[0])
    elif face_type == 'vi':
        path1a = (path_ids[1], path_ids[2])
        path1b = (path_ids[3], path_ids[0])
        path2a = (path_ids[2], path_ids[3])
        path2b = (path_ids[0], path_ids[1])
    elif face_type == 'vii':
        path1a = (path_ids[2], path_ids[3])
        path1b = (path_ids[0], path_ids[1])
        path2a = (path_ids[0], path_ids[3])
        path2b = (path_ids[2], path_ids[1])

    p1a_d = get_atom_distance(face, path1a[0], path1a[1])
    p1b_d = get_atom_distance(face, path1b[0], path1b[1])
    p2a_d = get_atom_distance(face, path2a[0], path2a[1])
    p2b_d = get_atom_distance(face, path2b[0], path2b[1])
    mismatch1 = (abs(p1a_d-p1b_d)/max([p1a_d, p1b_d])) * 100
    mismatch2 = (abs(p2a_d-p2b_d)/max([p2a_d, p2b_d])) * 100
    difference1 = abs(p1a_d-p1b_d)
    difference2 = abs(p2a_d-p2b_d)

    xys = (
        (
            (
                # X coordinates of vector 1.
                tuple(face.get_atomic_positions(path1a[0]))[0][0],
                tuple(face.get_atomic_positions(path1a[1]))[0][0],

            ),
            (
                # Y coordinates of vector 1.
                tuple(face.get_atomic_positions(path1a[0]))[0][1],
                tuple(face.get_atomic_positions(path1a[1]))[0][1],

            ),
            'r', '-',
        ),
        (
            (
                # X coordinates of vector 2.
                tuple(face.get_atomic_positions(path1b[0]))[0][0],
                tuple(face.get_atomic_positions(path1b[1]))[0][0],
            ),
            (
                # Y coordinates of vector 2.
                tuple(face.get_atomic_positions(path1b[0]))[0][1],
                tuple(face.get_atomic_positions(path1b[1]))[0][1],
            ),
            'skyblue', '-',
        ),
        (
            (
                # X coordinates of vector 1.
                tuple(face.get_atomic_positions(path2a[0]))[0][0],
                tuple(face.get_atomic_positions(path2a[1]))[0][0],

            ),
            (
                # Y coordinates of vector 1.
                tuple(face.get_atomic_positions(path2a[0]))[0][1],
                tuple(face.get_atomic_positions(path2a[1]))[0][1],

            ),
            'orange', '--',
        ),
        (
            (
                # X coordinates of vector 2.
                tuple(face.get_atomic_positions(path2b[0]))[0][0],
                tuple(face.get_atomic_positions(path2b[1]))[0][0],
            ),
            (
                # Y coordinates of vector 2.
                tuple(face.get_atomic_positions(path2b[0]))[0][1],
                tuple(face.get_atomic_positions(path2b[1]))[0][1],
            ),
            'green', '--',
        ),
    )

    return {
        'path1a': path1a,
        'path1b': path1b,
        'path2a': path2a,
        'path2b': path2b,
        'p1a_d': p1a_d,
        'p1b_d': p1b_d,
        'p2a_d': p2a_d,
        'p2b_d': p2b_d,
        'mismatch1': mismatch1,
        'mismatch2': mismatch2,
        'difference1': difference1,
        'difference2': difference2,
        'xys': xys,
    }


def visualise_face(face, face_name, face_type, paths):
    """
    Plot a visualisation of the face.

    """

    fig, ax = plt.subplots(figsize=(5, 5))
    plot_properties = {
        'metals': {'c': 'orange', 's': 160},
        'Cs': {'c': 'gray', 's': 120},
        'Ns': {'c': 'skyblue', 's': 120},
    }
    string = ''
    # Plot the properties of each path.
    xmin = 100
    xmax = -100
    ymin = 100
    ymax = -100
    for path in paths:
        if path != 'metals':
            continue
        path_atom_ids = paths[path]
        # Get mismatches.
        mismatch_values = calculate_path_mismatch(
            face=face,
            path_ids=path_atom_ids,
            face_type=face_type,
        )

        string += (
            f"{path}: ({round(mismatch_values['p1a_d'], 2)}, "
            f"{round(mismatch_values['p1b_d'], 2)}), "
            f"({round(mismatch_values['p2a_d'], 2)}, "
            f"{round(mismatch_values['p2b_d'], 2)}) "
            f"AR: {round(mismatch_values['mismatch1'], 2)}%, "
            f"{round(mismatch_values['mismatch2'], 2)}% "
            f"AB: {round(mismatch_values['difference1'], 2)}A, "
            f"{round(mismatch_values['difference2'], 2)}A\n"
        )

        # Plot paths.
        for xys in mismatch_values['xys']:
            ax.plot(
                xys[0], xys[1],
                c=xys[2],
                lw=2,
                linestyle=xys[3],
            )
        # Plot atom positions.
        for i in face.get_atomic_positions(path_atom_ids):
            ax.scatter(
                i[0], i[1],
                c=plot_properties[path]['c'],
                s=plot_properties[path]['s'],
                alpha=1.0,
                edgecolor='k',
            )
            xmax = max([xmax, i[0]])
            xmin = min([xmin, i[0]])
            ymax = max([ymax, i[1]])
            ymin = min([ymin, i[1]])

    xmid = 0  # (xmax+xmin)/2
    ymid = 0  # (ymax+ymin)/2
    ax.text(xmid-5, ymid-3, string, fontsize=8)

    # ax.tick_params(axis='both', which='major', labelsize=16)
    # ax.set_xlabel(r'$x$ [$\mathrm{\AA}}$]', fontsize=16)
    # ax.set_ylabel(r'$y$ [$\mathrm{\AA}}$]', fontsize=16)
    ax.set_aspect(1)
    ax.axis('off')
    filename = f'{face_name}_viz.pdf'
    fig.tight_layout()
    fig.savefig(filename, dpi=720, bbox_inches='tight')
    plt.close()


def get_face_properties(face, face_name, face_type, paths):

    json_file = f'{face_name}_properties.json'
    if exists(json_file):
        with open(json_file, 'r') as f:
            data = json.load(f)
    else:
        data = calculate_face_properties(face, paths, face_type)
        with open(json_file, 'w') as f:
            json.dump(data, f)

    return data


def face_convert(string):
    conv = {
        'i': 1,
        'ii': 2,
        'iii': 3,
        'iv': 4,
        'v': 5,
        'vi': 6,
        'vii': 7,
    }
    return conv[string]


def plot_face_mismatches(data, name, types='metals'):

    m1 = {}
    m2 = {}
    diff = {}
    avg = {}

    for face in data:
        m1[face] = data[face][types]['mms'][0]
        m2[face] = data[face][types]['mms'][1]
        # avg[face] = np.average(data[face][types]['mms'])
        # diff[face] = abs(
        #     data[face][types]['mms'][0] - data[face][types]['mms'][1]
        # )

    x_ticks = [face_convert(i) for i in data]
    x_ticklabels = [f'${i}$' for i in data]

    # width = 0.9
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.scatter(
        x=[face_convert(i)-0.1 for i in m1],
        y=[m1[i] for i in m1],
        # width=width,
        facecolor='gold',
        edgecolor='k',
        # linewidth=2,
        s=80,
        marker='o',
        alpha=1.0,
        # label='corner 1',
        # label=r'$M_\mathrm{F}$',
    )

    ax.scatter(
        x=[face_convert(i)+0.1 for i in m2],
        y=[m2[i] for i in m2],
        # width=width,
        facecolor='gold',
        edgecolor='k',
        # linewidth=2,
        s=80,
        marker='o',
        alpha=1.0,
        label=r'$M_\mathrm{F}$',
    )

    # ax.plot(
    #     [face_convert(i) for i in avg],
    #     [avg[i] for i in avg],
    #     # width=width,
    #     color='gold',
    #     # edgecolor='k',
    #     marker='o',
    #     markersize=9,
    #     linewidth=4,
    #     alpha=1,
    #     label='average',
    # )

    # ax.scatter(
    #     x=[face_convert(i) for i in diff],
    #     y=[diff[i] for i in diff],
    #     # width=width,
    #     facecolor='k',
    #     edgecolor='k',
    #     # linewidth=2,
    #     s=100,
    #     marker='D',
    #     alpha=1,
    #     label='difference',
    # )

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylabel('mismatch [%]', fontsize=16)
    ax.set_ylim(0, 60)
    # Set number of ticks for x-axis
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_ticklabels)
    ax.legend(fontsize=16)

    fig.tight_layout()
    if types == 'metals':
        fig.savefig(
            f'f_mismatch_{name}.pdf',
            dpi=720,
            bbox_inches='tight'
        )
    else:
        fig.savefig(
            f'f_mismatch_{types}_{name}.pdf',
            dpi=720,
            bbox_inches='tight'
        )
    plt.close()


def plot_face_differences(data, name, types='metals'):

    m1 = {}
    m2 = {}
    diff = {}
    avg = {}

    for face in data:
        m1[face] = data[face][types]['dif'][0]
        m2[face] = data[face][types]['dif'][1]
        # avg[face] = np.average(data[face][types]['dif'])
        # diff[face] = abs(
        #     data[face][types]['dif'][0] - data[face][types]['dif'][1]
        # )

    x_ticks = [face_convert(i) for i in data]
    x_ticklabels = [f'${i}$' for i in data]

    # width = 0.9
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.scatter(
        x=[face_convert(i)-0.1 for i in m1],
        y=[m1[i] for i in m1],
        # width=width,
        facecolor='gold',
        edgecolor='k',
        # linewidth=2,
        s=80,
        marker='o',
        alpha=1.0,
        # label='corner 1',
    )

    ax.scatter(
        x=[face_convert(i)+0.1 for i in m2],
        y=[m2[i] for i in m2],
        # width=width,
        facecolor='gold',
        edgecolor='k',
        # linewidth=2,
        s=80,
        marker='o',
        alpha=1.0,
        label=r'$D_\mathrm{F}$',
    )

    # ax.plot(
    #     [face_convert(i) for i in avg],
    #     [avg[i] for i in avg],
    #     # width=width,
    #     color='gold',
    #     # edgecolor='k',
    #     marker='o',
    #     markersize=9,
    #     linewidth=4,
    #     alpha=1,
    #     label='average',
    # )

    # ax.scatter(
    #     x=[face_convert(i) for i in diff],
    #     y=[diff[i] for i in diff],
    #     # width=width,
    #     facecolor='k',
    #     edgecolor='k',
    #     # linewidth=2,
    #     s=100,
    #     marker='D',
    #     alpha=1,
    #     label='difference',
    # )

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylabel('difference [$\mathrm{\AA}$]', fontsize=16)
    ax.set_ylim(0, 20)
    # Set number of ticks for x-axis
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_ticklabels)
    ax.legend(fontsize=16)

    fig.tight_layout()
    if types == 'metals':
        fig.savefig(
            f'f_differences_{name}.pdf',
            dpi=720,
            bbox_inches='tight'
        )
    else:
        fig.savefig(
            f'f_differences_{types}_{name}.pdf',
            dpi=720,
            bbox_inches='tight'
        )
    plt.close()


def get_paths(face, face_name, metal_atomic_number=30):
    """
    Get the atom id paths to use to calculate the mismatch.

    """

    # bb id: iter of neighbouring bb ids defined in CubeFace.
    # Match these with the expected coordinates for specific corners.
    bb_id_neigh_order = [None, None, None, None]
    # Also get metal atom ids.
    metal_atom_ids = []
    for atomi in face.get_atom_infos():
        atom = atomi.get_atom()
        if atom.get_atomic_number() == metal_atomic_number:
            metal_atom_ids.append(atom.get_id())
            bbid = atomi.get_building_block_id()
            pos = tuple(face.get_atomic_positions(atom.get_id()))[0]
            x, y, _ = pos
            if x > 0:
                if y > 0:
                    bb_id_neigh_order[0] = bbid
                else:
                    bb_id_neigh_order[1] = bbid
            else:
                if y > 0:
                    bb_id_neigh_order[3] = bbid
                else:
                    bb_id_neigh_order[2] = bbid

    if None in bb_id_neigh_order:
        print(bb_id_neigh_order)
        raise ValueError(
            f'BB ID neighbour search failed for {face_name}'
        )

    N_atom_ids = []
    C_atom_ids = []
    # N atom id: C atom id
    NC_pairs = {}
    smarts = '[#6X3]-[#7]=[#6X3H1]-[#6X3!H1]'
    rdkit_mol = face.to_rdkit_mol()
    query_ids = get_query_atom_ids(smarts, rdkit_mol)
    for atom_ids in query_ids:
        C_atom_ids.append(atom_ids[0])
        N_atom_ids.append(atom_ids[1])
        NC_pairs[atom_ids[1]] = atom_ids[0]
    if len(N_atom_ids) != 4:
        raise ValueError(f'too many matches found in {face_name}!')

    # Atom ids that make the path around the face matching the bb id
    # neighbour order.
    paths = {
        'metals': [0, 0, 0, 0],
        'Cs': [0, 0, 0, 0],
        'Ns': [0, 0, 0, 0],
    }
    for ai in face.get_atom_infos():
        atom = ai.get_atom()
        atom_id = atom.get_id()
        bb_id = ai.get_building_block_id()
        if atom_id in metal_atom_ids:
            order_idx = bb_id_neigh_order.index(bb_id)
            paths['metals'][order_idx] = atom_id

        if atom_id in N_atom_ids:
            order_idx = bb_id_neigh_order.index(bb_id)
            paths['Ns'][order_idx] = atom_id

    for i, N_id in enumerate(paths['Ns']):
        paired_C = NC_pairs[N_id]
        paths['Cs'][i] = paired_C

    return paths


def heatmap(
    data_dict,
    vmin,
    vmax,
):

    faces = ['i', 'ii', 'iii', 'iv', 'v']  # , 'vi', 'vii']
    _expt_lig_data = {
        'quad2_12': 'ii',
        'quad2_8': 'iii',
        'quad2_3': 'ii',
        'quad2_16': 'ii',
        'quad2_2': 'iii',
        'quad2_5': 'i',
    }

    fig, ax = plt.subplots(figsize=(8, 8))
    xshape = len(data_dict)
    yshape = len(faces)
    maps = np.zeros((xshape, yshape))
    for i, lig in enumerate(data_dict):
        da = data_dict[lig]
        print(da)
        raise SystemExit()
        for j, face in enumerate(da):
            maps[i][j] = np.average(da[face]['metals'])

    im = ax.imshow(maps, vmin=vmin, vmax=vmax, cmap='Purples_r')
    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, shrink=0.4)
    cbar.ax.set_ylabel(
        'avg. side mismatch [%]',
        rotation=-90, va="bottom", fontsize=16,
    )
    cbar.ax.tick_params(labelsize=16)

    # Min of each row.
    index_min = np.argmin(maps, axis=1)
    ax.scatter(
        x=index_min,
        y=[lig for lig in data_dict],
        c='white',
        marker='o',
        edgecolors='k',
        s=150,
    )

    for expt in _expt_lig_data:
        lig_position = list(data_dict.keys()).index(expt)
        face_position = faces.index(_expt_lig_data[expt])
        ax.scatter(
            x=face_position,
            y=lig_position,
            c='red',
            edgecolors='k',
            marker='P',
            s=120,
        )

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=1)
    # ax.set_xticks(np.arange(maps.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(maps.shape[0]+1)-.5, minor=True)
    ax.tick_params(which="minor", bottom=False, left=False)

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('face', fontsize=16)
    ax.set_ylabel('ligand', fontsize=16)

    # Show all ticks and label them with the respective lists.
    ax.set_xticks([face_convert(a)-1 for a in faces])
    ax.set_xticklabels([a for a in faces])
    ax.set_yticks([i for i in range(len(data_dict))])
    ax.set_yticklabels([lig for lig in data_dict])

    fig.tight_layout()
    fig.savefig(
        'face_map.pdf',
        dpi=720,
        bbox_inches='tight',
    )
    plt.close()


def main():
    first_line = (
        'Usage: face_analysis.py '
        'lig_directory complex_name complex_directory'
    )
    if (not len(sys.argv) == 4):
        print(f"""
{first_line}

    ligand_directory : (str)
        Directory with required ligand structures.

    complex_name : (str)
        Prefix of complex used in this run - defines face file names.

    complex_directory : (str)
        Location with manually constructed complexes.

    """)
        sys.exit()
    else:
        ligand_directory = sys.argv[1]
        complex_name = sys.argv[2]
        complex_directory = sys.argv[3]

    # Define two metal building blocks (lambda, delta).
    del_name, del_complex = load_complex(
        f'{complex_directory}/{complex_name}_zn_oct_del_face.mol'
    )
    lam_name, lam_complex = load_complex(
        f'{complex_directory}/{complex_name}_zn_oct_lam_face.mol'
    )

    # Optimise both complexes.
    del_complex = optimize_complex(del_complex, del_name)
    lam_complex = optimize_complex(lam_complex, lam_name)

    # Load in each ligand structure.
    ligands = load_ligands(ligand_directory)

    face_topologies = face_topology_dict()

    # Build and optimise five face options per ligand.
    face_matches = {}
    long_face_matches = {}
    for lig in sorted(ligands):
        print(f'doing {lig}...')
        lig_structure = ligands[lig]
        lig_structure.show_long_axis(
            long_axis=lig_structure.get_long_axis(),
            path=f'{lig}_lashown.xyz'
        )
        # Build each face topology.
        lig_faces = {}
        long_lig_faces = {}
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
            show_long_axis(face, face_name)
            opt_face = optimize_face(face, face_name)
            all_bls = get_all_bond_lengths(opt_face)
            if max(all_bls) > 3:
                raise ValueError(
                    f'max bond length of {face_name} is {max(all_bls)}'
                )

            # Get the paths to use in visualisation and calculation.
            paths = get_paths(opt_face, face_name)

            # Measure properties.
            face_properties = get_face_properties(
                opt_face, face_name, face_t, paths
            )
            print(
                f":: {face_name}: "
                f"{np.average(face_properties['metals']['dif'])}, "
                f"{np.average(face_properties['Ns']['dif'])}, "
                f"{np.average(face_properties['Cs']['dif'])}, "
            )
            lig_faces[face_t] = face_properties

            long_face_name = f'{face_name}_long'
            long_opt_face = long_optimize_face(
                face=opt_face,
                face_name=long_face_name,
            )
            all_bls = get_all_bond_lengths(long_opt_face)
            if max(all_bls) > 3:
                raise ValueError(
                    f'max bond length of {long_face_name} is '
                    f'{max(all_bls)}'
                )

            # Get the paths to use in visualisation and calculation.
            paths = get_paths(long_opt_face, long_face_name)
            # For visualisation.
            visualise_face(long_opt_face, face_name, face_t, paths)
            # Measure properties.
            long_face_properties = get_face_properties(
                long_opt_face, long_face_name, face_t, paths
            )
            print(
                f":: {long_face_name}: "
                f"{np.average(long_face_properties['metals']['dif'])}, "
                f"{np.average(long_face_properties['Ns']['dif'])}, "
                f"{np.average(long_face_properties['Cs']['dif'])}, "
            )
            long_lig_faces[face_t] = long_face_properties

        plot_face_mismatches(data=lig_faces, name=lig)
        face_matches[lig] = lig_faces
        plot_face_mismatches(data=long_lig_faces, name=f'{lig}_long')
        plot_face_differences(data=lig_faces, name=lig)
        plot_face_differences(data=long_lig_faces, name=f'{lig}_long')
        long_face_matches[lig] = long_lig_faces

    heatmap(
        data_dict=face_matches,
        vmin=0,
        vmax=40,
    )

    heatmap(
        data_dict=long_face_matches,
        vmin=0,
        vmax=40,
    )


if __name__ == '__main__':
    main()
