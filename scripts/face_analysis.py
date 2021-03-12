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

from atools import get_query_atom_ids

from molecule_building import metal_FFs
from cubeface import CubeFace, FaceBuildingBlock

from atools import (
    MOC_collapse_mc,
    MOC_uff_opt,
    get_atom_distance,
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
        ligands[l_name] = FaceBuildingBlock.init_from_file(
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


def calculate_face_properties(face, paths, metal_atomic_number=30):
    """
    Calculate geometrical properties of a face.

    """

    properties = {
        'metals': None,
        'Ns': None,
        'Cs': None,
    }
    for path in paths:
        path_atom_ids = paths[path]
        # Get mismatches.
        path1a = (path_atom_ids[1], path_atom_ids[0])
        path1b = (path_atom_ids[-1], path_atom_ids[0])
        path2a = (path_atom_ids[1], path_atom_ids[2])
        path2b = (path_atom_ids[-1], path_atom_ids[2])
        p1a_d = get_atom_distance(face, path1a[0], path1a[1])
        p1b_d = get_atom_distance(face, path1b[0], path1b[1])
        p2a_d = get_atom_distance(face, path2a[0], path2a[1])
        p2b_d = get_atom_distance(face, path2b[0], path2b[1])
        mismatch1 = (abs(p1a_d-p1b_d)/max([p1a_d, p1b_d])) * 100
        mismatch2 = (abs(p2a_d-p2b_d)/max([p2a_d, p2b_d])) * 100

        properties[path] = (mismatch1, mismatch2)

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


def visualise_face(face, face_name, paths):
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
        path_atom_ids = paths[path]
        # Get mismatches.
        path1a = (path_atom_ids[1], path_atom_ids[0])
        path1b = (path_atom_ids[-1], path_atom_ids[0])
        path2a = (path_atom_ids[1], path_atom_ids[2])
        path2b = (path_atom_ids[-1], path_atom_ids[2])
        p1a_d = get_atom_distance(face, path1a[0], path1a[1])
        p1b_d = get_atom_distance(face, path1b[0], path1b[1])
        p2a_d = get_atom_distance(face, path2a[0], path2a[1])
        p2b_d = get_atom_distance(face, path2b[0], path2b[1])
        mismatch1 = (abs(p1a_d-p1b_d)/max([p1a_d, p1b_d])) * 100
        mismatch2 = (abs(p2a_d-p2b_d)/max([p2a_d, p2b_d])) * 100

        string += (
            f'{path}: ({round(p1a_d, 2)}, {round(p1b_d, 2)}), '
            f'({round(p2a_d, 2)}, {round(p2b_d, 2)}) '
            f'AR: {round(mismatch1, 2)}%, {round(mismatch2, 2)}%\n'
        )

        # Plot paths.
        for pp in [path1a, path1b, path2a, path2b]:
            pos1 = tuple(face.get_atomic_positions(pp[0]))[0]
            pos2 = tuple(face.get_atomic_positions(pp[1]))[0]
            x1 = pos1[0]
            y1 = pos1[1]
            x2 = pos2[0]
            y2 = pos2[1]
            ax.plot(
                [x1, x2], [y1, y2],
                c=plot_properties[path]['c'],
                lw=2
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


def get_face_properties(face, face_name, paths):

    json_file = f'{face_name}_properties.json'
    if exists(json_file):
        with open(json_file, 'r') as f:
            data = json.load(f)
    else:
        data = calculate_face_properties(face, paths)
        with open(json_file, 'w') as f:
            json.dump(data, f)

    return data


def plot_face_mismatches(data, name):

    avg_m_mismatches = {}
    avg_n_mismatches = {}
    avg_c_mismatches = {}

    for face in data:
        avg_m_mismatch = np.average(data[face]['metals'])
        avg_n_mismatch = np.average(data[face]['Ns'])
        avg_c_mismatch = np.average(data[face]['Cs'])
        avg_m_mismatches[face] = avg_m_mismatch
        avg_n_mismatches[face] = avg_n_mismatch
        avg_c_mismatches[face] = avg_c_mismatch

    x_ticks = [int(i) for i in data]
    x_ticklabels = [i for i in data]

    # width = 0.9
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.scatter(
        x=[int(i) for i in avg_m_mismatches],
        y=[avg_m_mismatches[i] for i in avg_m_mismatches],
        # width=width,
        facecolor='orange',
        edgecolor='k',
        # linewidth=2,
        s=160,
        marker='o',
        alpha=1,
        label='M-M',
    )

    ax.scatter(
        x=[int(i) for i in avg_n_mismatches],
        y=[avg_n_mismatches[i] for i in avg_n_mismatches],
        # width=width,
        facecolor='skyblue',
        edgecolor='k',
        # linewidth=2,
        s=160,
        marker='P',
        alpha=1,
        label='N-N',
    )

    # ax.bar(
    #     x=[int(i) for i in avg_c_mismatches],
    #     height=[avg_c_mismatches[i] for i in avg_c_mismatches],
    #     width=width,
    #     facecolor='none',
    #     edgecolor='gray',
    #     linewidth=2,
    #     alpha=1,
    #     label='C-C',
    # )

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylabel('avg. face-side mismatch [%]', fontsize=16)
    ax.set_ylim(0, 75)
    # Set number of ticks for x-axis
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_ticklabels)
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(
        f'f_mismatch_{name}.pdf',
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
            'd_pos': (1, 3),
            'l_pos': (0, 2),
        },
        '2': {
            'va': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
            'd_pos': (0, 2),
            'l_pos': (1, 3),
        },
        '3': {
            'va': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
            'd_pos': (),
            'l_pos': (0, 1, 2, 3),
        },
        '4': {
            'va': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
            'd_pos': (3, ),
            'l_pos': (0, 1, 2),
        },
        '5': {
            'va': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
            'd_pos': (2, ),
            'l_pos': (0, 1, 3, ),
        },
        '6': {
            'va': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
            'd_pos': (0, 3),
            'l_pos': (1, 2),
        },
        '7': {
            'va': {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
            'd_pos': (0, 1),
            'l_pos': (2, 3),
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
        lig_structure.show_long_axis(
            long_axis=lig_structure.get_long_axis(),
            path=f'{lig}_lashown.xyz'
        )
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
            show_long_axis(face, face_name)
            opt_face = optimize_face(face, face_name)
            all_bls = get_all_bond_lengths(opt_face)
            if max(all_bls) > 3:
                raise ValueError(
                    f'max bond length of {face_name} is {max(all_bls)}'
                )

            # Get the paths to use in visualisation and calculation.
            paths = get_paths(opt_face, face_name)

            # For visualisation.
            visualise_face(opt_face, face_name, paths)

            # Measure properties.
            face_properties = get_face_properties(
                opt_face, face_name, paths
            )
            print(
                f":: {face_name}: "
                f"{np.average(face_properties['metals'])}, "
                f"{np.average(face_properties['Ns'])}, "
                f"{np.average(face_properties['Cs'])}, "
            )
            lig_faces[face_t] = face_properties
        plot_face_mismatches(data=lig_faces, name=lig)


if __name__ == '__main__':
    main()
