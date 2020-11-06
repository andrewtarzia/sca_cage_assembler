#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Utilities module.

Author: Andrew Tarzia

Date Created: 15 Mar 2020
"""

import numpy as np
from os.path import exists
import json
from itertools import combinations, permutations
import matplotlib.pyplot as plt
import matplotlib

from atools import (
    get_atom_distance,
    build_conformers,
    calculate_molecule_planarity,
    update_from_rdkit_conf,
    angle_between,
)


def read_lib(lib_file):
    """
    Read lib file.

    Returns dictionary.

    """

    print(f'reading {lib_file}')
    with open(lib_file, 'rb') as f:
        lib = json.load(f)

    return lib


def calculate_binding_AR(mol):
    """
    Calculate ligand aspect ratio based on binder positions.

    Defined as:
        Average ratio of the two shortest binding atom-binding atom
        distances eminating from each binding atom in the molecule.

    """

    if mol.get_num_functional_groups() != 4:
        return None

    binder_atom_ids = [
        list(fg.get_bonder_ids())
        for fg in mol.get_functional_groups()
    ]
    binder_atom_dists = sorted(
        [
            (idx1, idx2, get_atom_distance(mol, idx1, idx2))
            for idx1, idx2 in combinations(binder_atom_ids, r=2)
        ],
        key=lambda a: a[2]
    )

    far_binder_pair = (
        binder_atom_dists[-1][0],
        binder_atom_dists[-1][1]
    )
    ARs = []
    for fg_id in far_binder_pair:
        ds = sorted([
            i[2] for i in binder_atom_dists
            if fg_id in (i[0], i[1])
        ])
        AR = ds[1]/min(ds)
        ARs.append(AR)

    ligand_AR = sum(ARs)/len(ARs)
    return ligand_AR


def calculate_paired_face_anisotropies(mol, metal_atom_ids, face_sets):
    """
    Calculate face anisotropy of opposing sides of a prism.

    Defined as:
        Ratio of metal-metal distances defined by two edges eminating
        from a chosen metal on each face.

    """

    # Get all metal-metal distances with metal atom ids.
    metal_atom_dists = sorted(
        [
            (idx1, idx2, get_atom_distance(mol, idx1, idx2))
            for idx1, idx2 in combinations(metal_atom_ids, r=2)
        ],
        key=lambda a: a[2]
    )

    face_anisotropies = {}
    for fs in face_sets.sets:
        fsv = face_sets.vertices[fs]
        fsc = face_sets.connected[fs]
        fs_atom_ids = tuple(metal_atom_ids[i] for i in fsv)
        fs_distances = tuple(
            i for i in metal_atom_dists
            if i[0] in fs_atom_ids and i[1] in fs_atom_ids
        )
        # Pick one metal.
        idx = fsv[0]
        conn = [
            j
            for i in fsc if idx in i
            for j in i if j != idx
        ]
        fs_idx = metal_atom_ids[idx]
        fs_conn = [metal_atom_ids[i] for i in conn]

        # Define aniso based on the distances between its two
        # connections.
        pair1 = (fs_idx, fs_conn[0])
        pair2 = (fs_idx, fs_conn[1])
        d1, = (
            i[2]
            for i in fs_distances
            if i[0] in pair1 and i[1] in pair1
        )
        d2, = (
            i[2]
            for i in fs_distances
            if i[0] in pair2 and i[1] in pair2
        )
        face_aniso = d2 / d1 if d2 > d1 else d1 / d2
        face_anisotropies[fs] = face_aniso

    # Pair the face anisotropies of opposing faces.
    paired_face_anisotropies = [
        (i, j, face_anisotropies[i], face_anisotropies[j])
        for i, j in combinations(face_anisotropies, r=2)
        if face_sets.opposite[i] == j
    ]

    return paired_face_anisotropies


def calculate_cube_likeness(mol, metal_atom_ids, face_sets):
    """
    Calculate cube likeness a prism.

    Defined as:
        XXXXXXX

    """

    pos_mat = mol.get_position_matrix()

    metal_metal_vectors = {
        (idx1, idx2): pos_mat[idx2]-pos_mat[idx1]
        for idx1, idx2 in permutations(metal_atom_ids, r=2)
    }

    cube_likeness = {fs: {} for fs in face_sets.sets}
    for fs in face_sets.sets:
        print(fs)
        fsv = face_sets.vertices[fs]
        fsc = face_sets.connected[fs]
        fs_atom_ids = tuple(metal_atom_ids[i] for i in fsv)
        print(fs_atom_ids)

        # Calculate the interior angles based on connected metals.
        interior_angles = {}
        for idx in fsv:
            conn = [
                j
                for i in fsc if idx in i
                for j in i if j != idx
            ]
            fs_idx = metal_atom_ids[idx]
            fs_conn = [metal_atom_ids[i] for i in conn]

            # Define angle based on the two connections.
            pair1 = (fs_idx, fs_conn[0])
            pair2 = (fs_idx, fs_conn[1])
            # print(idx, pair1, pair2)
            vector1 = metal_metal_vectors[(pair1)]
            vector2 = metal_metal_vectors[(pair2)]
            interior_angle = np.degrees(
                angle_between(vector1, vector2)
            )
            print(vector1, vector2, interior_angle)
            interior_angles[idx] = interior_angle
        print(interior_angles, sum(interior_angles.values()))

        metal_plane_deviation = calculate_molecule_planarity(
            mol=mol,
            plane_ids=fs_atom_ids,
            atom_ids=fs_atom_ids,
        )
        print('>>>>', fs, metal_plane_deviation)
        cube_likeness[fs]['interior_angles'] = interior_angles
        cube_likeness[fs]['metal_PD'] = metal_plane_deviation

    return cube_likeness


def convert_symm_names(symm_name):

    new_names = {
        'o1': r'O',
        'th1': r'T$_{h, 1}$',
        'th2': r'T$_{h, 2}$',
        't1': r'T',
        's61': r'S$_{6, 1}$',
        's62': r'S$_{6, 2}$',
        'd31': r'D$_{3, 1}$',
        'd32': r'D$_{3, 2}$',
        'c2v': r'C$_{2h}$',
        'c2h': r'C$_{2v}$',
    }

    return new_names[symm_name]


def heatmap(
    data,
    row_labels,
    col_labels,
    ax=None,
    cbar_kw={},
    cbarlabel="",
    **kwargs
):
    """
    Create a heatmap from a numpy array and two lists of labels.

    From: https://matplotlib.org/3.1.1/gallery/
    images_contours_and_fields/image_annotated_heatmap.html
    #sphx-glr-gallery-images-contours-and-fields-image-
    annotated-heatmap-py

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is
        plotted.  If not provided, use current axes or create a new
        one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.
        Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.

    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(
        top=True,
        bottom=False,
        labeltop=True,
        labelbottom=False
    )

    # Rotate the tick labels and set their alignment.
    plt.setp(
        ax.get_xticklabels(),
        rotation=-30,
        ha="right",
        rotation_mode="anchor"
    )

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)
    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    # ax.grid(which="minor", color="k", linestyle='-', linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(
    im,
    data=None,
    valfmt="{x:.2f}",
    textcolors=["black", "white"],
    threshold=None,
    na_points=None,
    **textkw
):
    """
    A function to annotate a heatmap.

    From: https://matplotlib.org/3.1.1/gallery/
    images_contours_and_fields/image_annotated_heatmap.html
    #sphx-glr-gallery-images-contours-and-fields-image-
    annotated-heatmap-py

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.
        Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should
        either use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used
        for values below a threshold, the second for those above.
        Optional.
    threshold
        Value in data units according to which the colors from
        textcolors are applied.  If None (the default) uses the middle
        of the colormap as separation.  Optional.
    na_points
        Tuple of (row, column) indices whose text should be manually
        set to N/A. Optional
    **kwargs
        All other arguments are forwarded to each call to `text` used
        to create the text labels.

    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(
        horizontalalignment="center",
        verticalalignment="center"
    )
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(
                color=textcolors[int(im.norm(data[i, j]) > threshold)]
            )
            cond1 = na_points is not None
            cond2 = False if na_points is None else (i, j) in na_points
            if cond1 and cond2:
                text = im.axes.text(j, i, 'N/A', **kw)
            else:
                text = im.axes.text(
                    j,
                    i,
                    valfmt(data[i, j], None),
                    **kw
                )
            texts.append(text)

    return texts


def get_planar_conformer(molecule):
    cids, confs = build_conformers(
        mol=molecule,
        N=100,
        ETKDG_version='v3'
    )
    print(f'getting optimal conformer...')
    min_plane_dev = 100000000
    min_cid = -10

    new_molecule = molecule.clone()

    for cid in cids:

        # Update stk_mol to conformer geometry.
        new_molecule = update_from_rdkit_conf(
            stk_mol=new_molecule,
            rdk_mol=confs,
            conf_id=cid
        )

        plane_dev = calculate_molecule_planarity(new_molecule)
        if plane_dev < min_plane_dev:
            min_cid = cid
            min_plane_dev = plane_dev
            molecule = update_from_rdkit_conf(
                stk_mol=molecule,
                rdk_mol=confs,
                conf_id=min_cid
            )

    return molecule


def planarfy(ligands):
    """
    Get the most planar conformer of each ligand.

    This is done by determining the ETKDG conformer with the smallest
    plane deviation from its plane of best fit.

    """

    new_ligands = {}

    for ligand in ligands:
        planar_file = f'{ligand}_planar.mol'
        if exists(planar_file):
            opt_lig = ligands[ligand].with_structure_from_file(
                planar_file
            )
        else:
            print(f'doing {ligand}...')
            opt_lig = get_planar_conformer(ligands[ligand])
            opt_lig.write(planar_file)
        new_ligands[ligand] = opt_lig

    return new_ligands


def split_xyz_file(num_atoms, xyz_file):
    """
    Splits xyz trajectory file into xyz files.

    """

    with open(xyz_file, 'r') as f:
        lines = f.readlines()

    file_strings = []
    string = []
    for line in lines:
        if f' {num_atoms} ' in f' {line.strip()} ':
            if len(string) == 0:
                string.append(line)
            else:
                # New block.
                file_strings.append(string)
                string = [line]
        else:
            string.append(line)
    # Add last set.
    file_strings.append(string)

    out_files = []
    for i, fs in enumerate(file_strings):
        file_name = xyz_file.replace('.xyz', f'_s{i}.xyz')
        with open(file_name, 'w') as f:
            for line in fs:
                f.write(line)
        out_files.append(file_name)

    return out_files
