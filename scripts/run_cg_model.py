#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run CG model analysis.

Author: Andrew Tarzia

Date Created: 17 Feb 2022

"""

import sys
import stk
import numpy as np
import rdkit.Chem.AllChem as rdkit
import os
import re
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import json
from string import digits

import stk

from env_set import gulp_path, shape_path
from symmetries import M8L6_Symmetry
from utilities import (
    reorient_linker,
    run_shape,
    ref_shape_dict,
    collect_all_shape_values,
    convert_symm_names,
    read_lib,
    angle_between,
    get_atom_distance,
)
from facebuildingblock import FaceBuildingBlock


def get_all_angles(molecule):

    paths = rdkit.FindAllPathsOfLengthN(
        mol=molecule.to_rdkit_mol(),
        length=3,
        useBonds=False,
        useHs=True,
    )
    angles = []
    for atom_ids in paths:
        atoms = list(
            molecule.get_atoms(atom_ids=[i for i in atom_ids])
        )
        atom1 = atoms[0]
        atom2 = atoms[1]
        atom3 = atoms[2]
        angles.append((atom1, atom2, atom3))

    return angles


class CGGulpOptimizer:

    def __init__(
        self,
        fileprefix,
        output_dir,
        anisotropy,
        ortho_k,
        o_angle_k,
    ):
        self._fileprefix = fileprefix
        self._output_dir = output_dir
        self._anisotropy = anisotropy
        self._ortho_k = ortho_k
        self._o_angle_k = o_angle_k
        self._gulp_in = os.path.join(
            self._output_dir, f'{self._fileprefix}.gin'
        )
        self._gulp_out = os.path.join(
            self._output_dir, f'{self._fileprefix}.ginout'
        )
        self._output_xyz = os.path.join(
            self._output_dir, f'{self._fileprefix}_final.xyz'
        )
        self._mass = 1
        self._bond_cutoff = 30
        self._angle_cutoff = 30

    def _run_gulp(self):
        os.system(
            f'{gulp_path()} < {self._gulp_in} > {self._gulp_out}'
        )

    def _extract_gulp(self):
        with open(self._gulp_out, 'r') as f:
            lines = f.readlines()

        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        run_data = {'traj': {}}
        for line in lines:
            if 'Cycle' in line:
                splits = line.rstrip().split()
                run_data['traj'][int(splits[1])] = {
                    'energy': float(splits[3]),
                    'gnorm': float(splits[5]),
                }

            if 'Final energy' in line:
                string = nums.search(line.rstrip()).group(0)
                energy = float(string)
                run_data['final_energy'] = energy

            if 'Final Gnorm' in line:
                string = nums.search(line.rstrip()).group(0)
                gnorm = float(string)
                run_data['final_gnorm'] = gnorm

        return run_data

    def define_bond_potentials(self):
        bond_ks_ = {
            ('C', 'C'): 10,
            ('B', 'B'): 10,
            ('B', 'C'): self._ortho_k,

            ('C', 'Zn'): 10,
            ('B', 'Zn'): 10,

            ('Fe', 'Fe'): 10,
            ('Fe', 'Zn'): 10,
            ('Zn', 'Zn'): 10,
        }
        _base_length = 4
        _ortho_length = self._anisotropy*_base_length
        bond_rs_ = {
            ('C', 'C'): _base_length,
            ('B', 'B'): _base_length,
            ('B', 'C'): _ortho_length,

            ('C', 'Zn'): 4,
            ('B', 'Zn'): 4,

            ('Fe', 'Fe'): 4,
            ('Fe', 'Zn'): 4,
            ('Zn', 'Zn'): 4,
        }
        return bond_ks_, bond_rs_

    def define_angle_potentials(self):
        angle_ks_ = {
            ('B', 'C', 'C'): self._o_angle_k,
            ('B', 'B', 'C'): self._o_angle_k,

            ('Fe', 'Fe', 'Fe'): 20,
            # ('Fe', 'Fe', 'Zn'): 10,
            ('Fe', 'Zn', 'Zn'): 20,
            ('Zn', 'Zn', 'Zn'): 20,

            ('B', 'C', 'Zn'): self._o_angle_k,
            ('B', 'B', 'Zn'): self._o_angle_k,
            ('C', 'C', 'Zn'): self._o_angle_k,

            ('C', 'Fe', 'Zn'): 20,
            ('B', 'Fe', 'Zn'): 20,
        }
        angle_thetas_ = {
            ('B', 'C', 'C'): 90,
            ('B', 'B', 'C'): 90,

            ('Fe', 'Fe', 'Fe'): 60,
            # ('Fe', 'Fe', 'Zn'): 60,
            # This requires a special rule.
            ('Fe', 'Zn', 'Zn'): (
                'check',
                {'cut': 70, 'min': 60, 'max': 90},
            ),
            ('Zn', 'Zn', 'Zn'): 60,

            ('B', 'C', 'Zn'): 135,
            ('B', 'B', 'Zn'): 135,
            ('C', 'C', 'Zn'): 135,

            ('C', 'Fe', 'Zn'): 180,
            ('B', 'Fe', 'Zn'): 180,
        }
        return angle_ks_, angle_thetas_

    def _get_coord_mass_string(self, mol):
        coord_string = 'cartesian\n'
        mass_string = ''

        pos_mat = mol.get_position_matrix()
        atoms = list(mol.get_atoms())
        for atom, pos_ in zip(atoms, pos_mat):
            name = f'{atom.__class__.__name__}{atom.get_id()+1}'
            coord_string += (
                f'{name} {round(pos_[0], 2)} {round(pos_[1], 2)} '
                f'{round(pos_[2], 2)}\n'
            )
            mass_string += f'mass {name} {self._mass}\n'

        return coord_string, mass_string

    def _get_bond_string(self, mol):
        bond_ks_, bond_rs_ = self.define_bond_potentials()
        bond_string = 'harm\n'
        bonds = list(mol.get_bonds())

        for bond in bonds:
            atom1 = bond.get_atom1()
            name1 = f'{atom1.__class__.__name__}{atom1.get_id()+1}'
            atom2 = bond.get_atom2()
            name2 = f'{atom2.__class__.__name__}{atom2.get_id()+1}'
            table = str.maketrans('', '', digits)
            sorted_name = tuple(sorted(
                [
                    i.translate(table)
                    for i in (name1, name2)
                ]
            ))

            try:
                bond_k = bond_ks_[sorted_name]
                bond_r = bond_rs_[sorted_name]
            except KeyError:
                continue

            bond_string += (
                f'{name1} {name2}  {bond_k} {bond_r} '
                f'{self._bond_cutoff}\n'
            )
        return bond_string

    def _get_angle_string(self, mol):
        angle_string = 'three\n'
        angle_ks_, angle_thetas_ = self.define_angle_potentials()
        angles = get_all_angles(mol)
        pos_mat = mol.get_position_matrix()

        for angle in angles:
            atom1, atom2, atom3 = angle
            name1 = f'{atom1.__class__.__name__}{atom1.get_id()+1}'
            name2 = f'{atom2.__class__.__name__}{atom2.get_id()+1}'
            name3 = f'{atom3.__class__.__name__}{atom3.get_id()+1}'
            table = str.maketrans('', '', digits)
            sorted_name = tuple(sorted(
                [
                    i.translate(table)
                    for i in (name1, name2, name3)
                ]
            ))

            try:
                angle_k = angle_ks_[sorted_name]
                angle_theta = angle_thetas_[sorted_name]
                if isinstance(angle_theta, int):
                    pass
                elif angle_theta[0] == 'check':
                    a1id = atom1.get_id()
                    a2id = atom2.get_id()
                    a3id = atom3.get_id()
                    vector1 = pos_mat[a2id]-pos_mat[a1id]
                    vector2 = pos_mat[a2id]-pos_mat[a3id]
                    curr_angle = np.degrees(
                        angle_between(vector1, vector2)
                    )
                    if curr_angle < angle_theta[1]['cut']:
                        angle_theta = angle_theta[1]['min']
                    elif curr_angle >= angle_theta[1]['cut']:
                        angle_theta = angle_theta[1]['max']

            except KeyError:
                continue

            angle_string += (
                f'{name2} {name1} {name3} {angle_k} {angle_theta} '
                f'{self._angle_cutoff} {self._angle_cutoff} '
                f'{self._angle_cutoff} \n'
            )

        return angle_string

    def _write_gulp_input(self, mol):
        top_string = 'opti conv cartesian\n'
        coord_string, mass_string = self._get_coord_mass_string(mol)
        bond_string = self._get_bond_string(mol)
        angle_string = self._get_angle_string(mol)
        settings_string = (
            '\nmaxcyc 500\n'
            # f'output xyz movie {filename}_traj.xyz\n'
            f'output xyz {self._output_xyz}\n'
        )

        with open(self._gulp_in, 'w') as f:
            f.write(top_string)
            f.write(coord_string)
            f.write(mass_string)
            f.write(bond_string)
            f.write(angle_string)
            f.write(settings_string)

    def optimize(self, molecule):
        self._write_gulp_input(mol=molecule)
        self._run_gulp()
        return self._extract_gulp()


def symmetries(cdelta, clambda, plane):
    symm_list = {}
    # Predefined list of symmetries.
    symm_c = M8L6_Symmetry(
        D_complex=cdelta,
        L_complex=clambda,
        linker=plane,
    )
    symm_list['d2'] = symm_c.d2()
    symm_list['th1'] = symm_c.th1()
    symm_list['th2'] = symm_c.th2()
    symm_list['td'] = symm_c.td()
    symm_list['tl'] = symm_c.tl()
    symm_list['s41'] = symm_c.s41()
    symm_list['s42'] = symm_c.s42()
    symm_list['s61'] = symm_c.s61()
    symm_list['s62'] = symm_c.s62()
    symm_list['d31'] = symm_c.d31()
    symm_list['d32'] = symm_c.d32()
    symm_list['d31n'] = symm_c.d31n()
    symm_list['d32n'] = symm_c.d32n()
    symm_list['c2v'] = symm_c.c2v()
    symm_list['c2h'] = symm_c.c2h()

    return symm_list


class CGM8L6Cube(stk.cage.M8L6Cube):

    def _get_scale(self, building_block_vertices):
        return 10


def write_cg_shape_input_file(
    input_file,
    structure_string,
    num_vertices,
    central_atom_id,
    ref_shapes,
):
    """
    Write input file for shape.

    """

    title = '$shape run by Andrew Tarzia.\n'
    size_of_poly = f'{num_vertices} {central_atom_id}\n'
    codes = ' '.join(ref_shapes)+'\n'

    string = title+size_of_poly+codes+structure_string

    with open(input_file, 'w') as f:
        f.write(string)


def calculate_cgcube_shape_measure(name, structure_string):
    """
    Calculate the shape of an 8 atom molecule.

    Shape: http://www.ee.ub.edu/index.php?option=com_content&view=
    article&id=575:shape-available&catid=80:news&Itemid=466

    """

    shape_dicts = (ref_shape_dict()['cube'], )
    n_verts = list(set([i['vertices'] for i in shape_dicts]))
    if len(n_verts) != 1:
        raise ValueError('Different vertex shapes selected.')

    input_file = f'{name}_shp.dat'
    std_out = f'{name}_shp.out'
    output_file = f'{name}_shp.tab'
    write_cg_shape_input_file(
        input_file=input_file,
        structure_string=structure_string,
        num_vertices=n_verts[0],
        central_atom_id=0,
        ref_shapes=[i['code'] for i in shape_dicts],
    )

    run_shape(input_file, shape_path(), std_out)
    shapes = collect_all_shape_values(output_file)
    return shapes


def prepare_precursors(precursor_dir):
    delta_bb = stk.BuildingBlock.init_from_file(
        path=os.path.join(precursor_dir, 'corner_delta.mol'),
        functional_groups=(stk.BromoFactory(), )
    )
    lambda_bb = stk.BuildingBlock.init_from_file(
        path=os.path.join(precursor_dir, 'corner_lambda.mol'),
        functional_groups=(stk.BromoFactory(), )
    )

    plane_bb = FaceBuildingBlock.init_from_file(
        path=os.path.join(precursor_dir, 'plane.mol'),
        functional_groups=(stk.BromoFactory(), )
    )

    temp_plane = reorient_linker(plane_bb)

    # Set functional group ordering based on long axis.
    fg_centroids = tuple(
        temp_plane.get_centroid(
            atom_ids=fg.get_placer_ids(),
        ) for fg in temp_plane.get_functional_groups()
    )
    plus_minus_fg_id = tuple(
        i for i, cent in enumerate(fg_centroids)
        if cent[0] > 0 and cent[1] < 0
    )[0]
    fg1_id = plus_minus_fg_id
    fg2_id, fg3_id, fg4_id = tuple(
        i
        for i in range(temp_plane.get_num_functional_groups())
        if i != fg1_id
    )
    new_fgs = tuple(temp_plane.get_functional_groups())
    plane_bb = temp_plane.with_functional_groups(
        functional_groups=(
            new_fgs[fg1_id],
            new_fgs[fg2_id],
            new_fgs[fg3_id],
            new_fgs[fg4_id],
        )
    )

    return delta_bb, lambda_bb, plane_bb


def get_shape_measure(cage, run_prefix, output_dir):
    Zn_bb_ids = {}
    for ai in cage.get_atom_infos():
        aibbid = ai.get_building_block_id()
        if ai.get_atom().get_atomic_number() == 30:
            if aibbid not in Zn_bb_ids:
                Zn_bb_ids[aibbid] = []
            Zn_bb_ids[aibbid].append(
                ai.get_atom().get_id()
            )

    Zn_centroids = []
    for n in Zn_bb_ids:
        Zn_centroids.append(cage.get_centroid(
            atom_ids=Zn_bb_ids[n]
        ))
    with open(
        os.path.join(output_dir, f'{run_prefix}_cents.xyz'), 'w'
    ) as f:
        f.write('8\n\n')
        for c in Zn_centroids:
            f.write(f'Zn {c[0]} {c[1]} {c[2]}\n')

    # Run calculations.
    s_string = f'{run_prefix}\n'
    for c in Zn_centroids:
        s_string += f'Zn {c[0]} {c[1]} {c[2]}\n'

    cu8_measure = calculate_cgcube_shape_measure(
        name=os.path.join(output_dir, f'{run_prefix}'),
        structure_string=s_string,
    )
    return cu8_measure


def get_distances(optimizer, cage):
    bond_ks_, __ = optimizer.define_bond_potentials()
    set_ks = tuple(bond_ks_.keys())
    distances = {''.join(i): [] for i in set_ks}
    for bond in cage.get_bonds():
        a1 = bond.get_atom1()
        a2 = bond.get_atom2()
        a1name = a1.__class__.__name__
        a2name = a2.__class__.__name__
        pair = tuple(sorted([a1name, a2name]))
        if pair in set_ks:
            a1id = a1.get_id()
            a2id = a2.get_id()
            distances[''.join(pair)].append(
                get_atom_distance(cage, a1id, a2id)
            )

    return distances


def get_angles(optimizer, cage):
    angle_ks_, __ = optimizer.define_angle_potentials()
    set_ks = tuple(angle_ks_.keys())
    angles = {''.join(i): [] for i in set_ks}
    pos_mat = cage.get_position_matrix()

    angle_atoms = get_all_angles(cage)
    for angle_trip in angle_atoms:
        triplet = tuple(
            sorted([i.__class__.__name__ for i in angle_trip])
        )
        if triplet in set_ks:
            a1id = angle_trip[0].get_id()
            a2id = angle_trip[1].get_id()
            a3id = angle_trip[2].get_id()
            vector1 = pos_mat[a2id]-pos_mat[a1id]
            vector2 = pos_mat[a2id]-pos_mat[a3id]
            angles[''.join(triplet)].append(np.degrees(
                angle_between(vector1, vector2)
            ))

    return angles


def run_aniso_optimisation(
    cage,
    aniso,
    symm,
    flex,
    ortho_k,
    o_angle_k,
    output_dir,
):

    run_prefix = f'{symm}_{aniso}_{flex}'
    output_file = os.path.join(
        output_dir, f'{run_prefix}_res.json'
    )

    if os.path.exists(output_file):
        with open(output_file, 'r') as f:
            res_dict = json.load(f)
    else:
        print(f': running optimisation of {run_prefix}')
        opt = CGGulpOptimizer(
            fileprefix=run_prefix,
            output_dir=output_dir,
            anisotropy=aniso,
            ortho_k=ortho_k,
            o_angle_k=o_angle_k,
        )
        run_data = opt.optimize(cage)

        # Get cube shape measure.
        opted = cage.with_structure_from_file(
            path=os.path.join(output_dir, f'{run_prefix}_final.xyz'),
        )
        opted.write(
            os.path.join(output_dir, f'{run_prefix}_final.mol')
        )

        cu8_measure = get_shape_measure(opted, run_prefix, output_dir)
        distances = get_distances(optimizer=opt, cage=opted)
        angles = get_angles(optimizer=opt, cage=opted)

        num_steps = len(run_data['traj'])
        fin_energy = run_data['final_energy']
        fin_gnorm = run_data['final_gnorm']
        traj_data = run_data['traj']
        print(
            f'{run_prefix}: {num_steps} {fin_energy} {fin_gnorm} '
            f'{cu8_measure}'
        )
        res_dict = {
            'fin_energy': fin_energy,
            'cu8': cu8_measure,
            'traj': traj_data,
            'distances': distances,
            'angles': angles,
        }
        with open(output_file, 'w') as f:
            json.dump(res_dict, f)

    return res_dict


def scatter(
    symm_to_c,
    results,
    ylabel,
    output_dir,
    filename,
    flex,
):

    fig, ax = plt.subplots(figsize=(8, 5))
    for aniso in results:
        da = results[aniso]
        for symm in da:
            if ylabel == 'energy (eV)':
                ys = da[symm]['fin_energy']
            elif ylabel == 'CU-8':
                ys = da[symm]['cu8']['CU-8']
                ax.axhline(y=0, lw=2, c='k')

            ax.scatter(
                aniso,
                ys,
                c=symm_to_c[symm][1],
                marker=symm_to_c[symm][0],
                edgecolor='k',
                s=80,
                alpha=0.5,
            )

    legend_elements = []
    for s in symm_to_c:
        legend_elements.append(
            Line2D(
                [0],
                [0],
                color='w',
                marker=symm_to_c[s][0],
                label=convert_symm_names(s),
                markerfacecolor=symm_to_c[s][1],
                markersize=10,
                markeredgecolor='k',
            )
        )

    ax.legend(handles=legend_elements, fontsize=16, ncol=3)

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('anisotropy', fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_title(f'flex: {flex}', fontsize=16)

    fig.tight_layout()
    fig.savefig(
        os.path.join(output_dir, filename),
        dpi=720,
        bbox_inches='tight',
    )
    plt.close()


def ey_vs_shape(
    results,
    output_dir,
    filename,
):

    _to_plot = {
        'd2': ('o', 'k'),
        'th2': ('X', 'r'),
        's62': ('D', 'gold'),
        'd32': ('o', 'skyblue'),
    }

    fig, ax = plt.subplots(figsize=(8, 5))
    for symm in _to_plot:
        x_vals = []
        y_vals = []
        for aniso in results:
            da = results[aniso]
            x_vals.append(da[symm]['cu8']['CU-8'])
            y_vals.append(da[symm]['fin_energy'])

        ax.scatter(
            x_vals,
            y_vals,
            c=_to_plot[symm][1],
            marker=_to_plot[symm][0],
            edgecolor='k',
            s=100,
            alpha=1.0,
            label=convert_symm_names(symm),
        )

    ax.legend(fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('CU-8', fontsize=16)
    ax.set_ylabel('energy (eV)', fontsize=16)
    ax.set_xlim(0, 2)

    fig.tight_layout()
    fig.savefig(
        os.path.join(output_dir, filename),
        dpi=720,
        bbox_inches='tight',
    )
    plt.close()


def comp_scatter(
    symm_to_c,
    symm_set,
    results,
    ylabel,
    output_dir,
    filename,
    flex,
    ylim,
):

    fig, ax = plt.subplots(figsize=(8, 5))
    for aniso in results:
        da = results[aniso]
        for symm in da:
            if symm not in symm_set:
                continue
            if ylabel == 'energy (eV)':
                ys = da[symm]['fin_energy']
            elif ylabel == 'CU-8':
                ys = da[symm]['cu8']['CU-8']
                ax.axhline(y=0, lw=2, c='k')

            ax.scatter(
                aniso,
                ys,
                c=symm_to_c[symm][1],
                marker=symm_to_c[symm][0],
                edgecolor='k',
                s=120,
            )

    legend_elements = []
    for s in symm_to_c:
        if s not in symm_set:
            continue
        legend_elements.append(
            Line2D(
                [0],
                [0],
                color='w',
                marker=symm_to_c[s][0],
                label=convert_symm_names(s),
                markerfacecolor=symm_to_c[s][1],
                markersize=12,
                markeredgecolor='k',
            )
        )

    ax.legend(handles=legend_elements, fontsize=16, ncol=2)

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('anisotropy', fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_title(f'flex: {flex}', fontsize=16)
    ax.set_ylim(ylim)

    fig.tight_layout()
    fig.savefig(
        os.path.join(output_dir, filename),
        dpi=720,
        bbox_inches='tight',
    )
    plt.close()


def merge_bond_types(s):

    translation = {
        'CC': 'face',
        'BB': 'face',
        'CZn': 'face-metal',
        'BZn': 'face-metal',
        'FeFe': 'metal',
        'FeZn': 'metal',
        'ZnZn': 'metal',
    }

    return translation[s]


def merge_angle_types(s):

    translation = {
        'BCC': 'face',
        'BBC': 'face',
        'BCZn': 'face-metal',
        'BBZn': 'face-metal',
        'CCZn': 'face-metal',
        'BFeZn': 'face-metal',
        'CFeZn': 'face-metal',
        'FeZnZn': 'metal',
        'FeFeFe': 'metal',
        'ZnZnZn': 'metal',
    }

    return translation[s]


def geom_distributions(
    results,
    output_dir,
    filename,
):

    # Collect all values for each bond and angle type.
    distance_by_type = {}
    angle_by_type = {}
    for aniso in results:
        da = results[aniso]
        for symm in da:
            dists = da[symm]['distances']
            angles = da[symm]['angles']
            for d in dists:
                if d == 'BC':
                    dd = f'{d}{aniso}'
                else:
                    dd = merge_bond_types(d)
                if dd in distance_by_type:
                    distance_by_type[dd].extend(dists[d])
                else:
                    distance_by_type[dd] = dists[d]
            for a in angles:
                aa = merge_angle_types(a)
                if aa in angle_by_type:
                    angle_by_type[aa].extend(angles[a])
                else:
                    angle_by_type[aa] = angles[a]

    fig, axs = plt.subplots(
        nrows=3,
        ncols=1,
        figsize=(8, 8),
    )
    # Plot distributions of each bond type.
    for btype in distance_by_type:
        if 'BC' in btype:
            continue
        data = distance_by_type[btype]
        axs[0].hist(
            x=data,
            bins=50,
            range=(3.6, 4.4),
            density=True,
            histtype='step',
            # color='',
            label=btype,
            lw=3,
        )
    axs[0].tick_params(axis='both', which='major', labelsize=16)
    axs[0].set_xlabel('distance [$\mathrm{\AA}}$]', fontsize=16)
    axs[0].set_ylabel('frequency', fontsize=16)
    axs[0].legend(fontsize=16, ncol=1)

    # Plot distributions of each variable bond type.
    for btype in distance_by_type:
        if 'BC' not in btype:
            continue
        data = distance_by_type[btype]
        aniso = float(btype.replace('BC', ''))
        axs[1].scatter(
            x=[aniso for i in data],
            y=data,
            color='gray',
            s=30,
            alpha=0.3,
            rasterized=True,
        )
    axs[1].tick_params(axis='both', which='major', labelsize=16)
    axs[1].set_xlabel('anisotropy', fontsize=16)
    axs[1].set_ylabel('F1-F2 distance [$\mathrm{\AA}}$]', fontsize=16)

    # Plot distributions of each angle type.
    for atype in angle_by_type:
        data = angle_by_type[atype]
        axs[2].hist(
            x=data,
            bins=50,
            range=(20, 182),
            density=True,
            histtype='step',
            # color='',
            label=atype,
            lw=3,
        )
    axs[2].tick_params(axis='both', which='major', labelsize=16)
    axs[2].set_xlabel('angle [degrees]', fontsize=16)
    axs[2].set_ylabel('frequency', fontsize=16)
    axs[2].legend(fontsize=16, ncol=1)

    fig.tight_layout()
    fig.savefig(
        os.path.join(output_dir, filename),
        dpi=720,
        bbox_inches='tight',
    )
    plt.close()


def heatmap(
    symm_to_c,
    results,
    output_dir,
    filename,
    vmin,
    vmax,
    clabel,
    flex,
    expt_data,
    ligand_ars,
):

    fig, ax = plt.subplots(figsize=(8, 8))
    maps = np.zeros((len(symm_to_c), len(results)))
    for j, aniso in enumerate(results):
        da = results[aniso]
        for i, symm in enumerate(symm_to_c):
            if clabel == 'energy (eV)':
                maps[i][j] = da[symm]['fin_energy']
            elif clabel == 'CU-8':
                maps[i][j] = da[symm]['cu8']['CU-8']

    im = ax.imshow(maps, vmin=vmin, vmax=vmax, cmap='Purples_r')
    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, shrink=0.4)
    cbar.ax.set_ylabel(clabel, rotation=-90, va="bottom", fontsize=16)
    cbar.ax.tick_params(labelsize=16)

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=1)
    # ax.set_xticks(np.arange(maps.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(maps.shape[0]+1)-.5, minor=True)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Scatter points for where experiments land.
    for expt in expt_data:
        da = expt_data[expt]
        # Get symm position.
        symm_position = symm_to_c[da['symmetry']][2]
        known_aniso = ligand_ars[da['ligand_name']]['N']
        sub_arr = [(a-known_aniso)**2 for a in results]
        matched_x_position = np.argmin(sub_arr)
        ax.scatter(
            x=matched_x_position,
            y=symm_position,
            c='red',
            edgecolors='k',
            marker='o',
            s=80,
        )

    # Scatter points for lowest/highest energy for each symm.
    if clabel == 'energy (eV)':
        # Min of each row.
        index_min = np.argmin(maps, axis=1)
        ax.scatter(
            x=index_min,
            y=[symm_to_c[symm][2] for symm in symm_to_c],
            c='white',
            marker='P',
            edgecolors='k',
            s=80,
        )
        # # Max of each row.
        # index_max = np.argmax(maps, axis=1)
        # ax.scatter(
        #     x=index_max,
        #     y=[symm_to_c[symm][2] for symm in symm_to_c],
        #     c='white',
        #     edgecolors='k',
        #     marker='X',
        #     s=40,
        # )

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('anisotropy', fontsize=16)
    ax.set_ylabel('symmetry', fontsize=16)
    # Show all ticks and label them with the respective lists.
    ax.set_xticks([i for i in range(len(results))])
    ax.set_xticklabels([a for a in results])
    ax.set_yticks([symm_to_c[symm][2] for symm in symm_to_c])
    ax.set_yticklabels([
        convert_symm_names(symm) for symm in symm_to_c
    ])

    # Rotate the tick labels and set their alignment.
    plt.setp(
        ax.get_xticklabels(),
        rotation=45,
        ha="right",
        rotation_mode="anchor",
    )
    ax.set_title(f'flex: {flex}', fontsize=16)

    fig.tight_layout()
    fig.savefig(
        os.path.join(output_dir, filename),
        dpi=720,
        bbox_inches='tight',
    )
    plt.close()


def convergence(
    results,
    output_dir,
    filename,
):

    # Pick examples to plot.
    _to_plot = (
        'td_2.0',
        # 'c2v_1.65',
        'c2v_1.5',
        # 'th2_1.6',
        'th2_1.1',
        # 's61_1.85',
        's61_1.3',
        's42_1.05',
    )

    fig, axs = plt.subplots(
        nrows=len(_to_plot),
        ncols=1,
        sharex=True,
        figsize=(8, 10),
    )

    for name, ax in zip(_to_plot, axs):
        symm, aniso = name.split('_')
        da = results[float(aniso)]
        traj = da[symm]['traj']
        traj_x = [i for i in traj]
        traj_e = [traj[i]['energy'] for i in traj]
        traj_g = [traj[i]['gnorm'] for i in traj]

        color = 'tab:red'
        ax.plot(
            traj_x,
            traj_e,
            lw=5,
            # marker='o',
            # markersize=12,
            color=color,
        )
        ax.tick_params(axis='y', labelcolor=color, labelsize=16)
        ax.set_yscale('log')

        # instantiate a second axes that shares the same x-axis
        ax2 = ax.twinx()
        color = 'tab:blue'
        ax2.plot(
            traj_x,
            traj_g,
            lw=5,
            # marker='X',
            # markersize=12,
            color=color,
        )
        ax2.tick_params(axis='y', labelcolor=color, labelsize=16)
        ax2.set_yscale('log')

        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.text(
            x=340, y=1100,
            s=f'{convert_symm_names(symm)}, aniso={aniso}',
            fontsize=16,
        )
        if name == 'th2_1.1':
            ax.set_ylabel('energy [eV]', fontsize=16)
            ax2.set_ylabel('Gnorm', fontsize=16)

    ax.set_xlabel('step', fontsize=16)
    ax.set_xticks(range(0, 501, 50))
    ax.set_xticklabels([str(i) for i in range(0, 501, 50)])

    fig.tight_layout()
    fig.savefig(
        os.path.join(output_dir, filename),
        dpi=720,
        bbox_inches='tight',
    )
    plt.close()


def get_ligand_ars(ligand_directory):
    json_file = os.path.join(ligand_directory, 'ligand_ARs.json')
    with open(json_file, 'r') as f:
        ar_data = json.load(f)
    return ar_data


def main():
    first_line = (
        'Usage: run_cg_model.py precursor_dir output_dir'
        ' expt_lib_file ligand_directory'
    )
    if (not len(sys.argv) == 5):
        print(f"""
{first_line}

    precursor_dir : (str)
        Directrory containing precursor structures.

    output_dir : (str)
        Directrory to output files to.

    expt_lib_file : (str)
        File containing experimental symmetry  information (XXXXX).

    ligand_directory : (str)
        Directory with required ligand structures.

    """)
        sys.exit()
    else:
        precursor_dir = sys.argv[1]
        output_dir = sys.argv[2]
        expt_lib_file = sys.argv[3]
        ligand_directory = sys.argv[4]

    expt_data = read_lib(expt_lib_file)

    # Get ligand aspect ratios.
    ligand_ars = get_ligand_ars(ligand_directory)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    delta_bb, lambda_bb, plane_bb = prepare_precursors(
        precursor_dir=precursor_dir,
    )

    # Make cage of each symmetry.
    symms = symmetries(delta_bb, lambda_bb, plane_bb)
    anisotropies = np.arange(1.0, 2.01, 0.05)
    # results = {i: {} for i in symms}
    flexes = {
        'low': (10, 20),
        'high': (0.1, 2.0),
    }
    for flex in flexes:
        results = {round(i, 2): {} for i in anisotropies}
        for symm in symms:
            topology_graph = CGM8L6Cube(
                building_blocks=symms[symm]['building_blocks'],
                vertex_alignments=symms[symm]['vertex_alignments'],
                num_processes=1,
            )
            cage = stk.ConstructedMolecule(topology_graph)
            cage.write(os.path.join(output_dir, f'{symm}_unopt.mol'))

            for aniso in anisotropies:
                aniso = round(aniso, 2)
                res_dict = run_aniso_optimisation(
                    cage=cage,
                    aniso=aniso,
                    symm=symm,
                    flex=flex,
                    ortho_k=flexes[flex][0],
                    o_angle_k=flexes[flex][1],
                    output_dir=output_dir,
                )
                results[aniso][symm] = res_dict

        symm_to_c = {
            'd2': ('o', 'k', 0),
            'th1': ('D', 'r', 1),
            'th2': ('X', 'r', 2),
            'td': ('o', 'r', 3),
            'tl': ('P', 'r', 4),
            's61': ('X', 'gold', 5),
            's62': ('D', 'gold', 6),
            's41': ('X', 'gray', 7),
            's42': ('D', 'gray', 8),
            'd31': ('P', 'skyblue', 9),
            'd32': ('o', 'skyblue', 10),
            'd31n': ('P', 'b', 11),
            'd32n': ('o', 'b', 12),
            'c2h': ('o', 'green', 13),
            'c2v': ('X', 'green', 14),
        }

        convergence(
            results=results,
            output_dir=output_dir,
            filename=f'convergence_{flex}.pdf',
        )

        ey_vs_shape(
            results=results,
            output_dir=output_dir,
            filename=f'e_vs_shape_{flex}.pdf',
        )

        geom_distributions(
            results=results,
            output_dir=output_dir,
            filename=f'dist_{flex}.pdf',
        )

        heatmap(
            symm_to_c=symm_to_c,
            results=results,
            output_dir=output_dir,
            filename=f'energy_map_{flex}.pdf',
            vmin=0,
            vmax=45,
            clabel='energy (eV)',
            flex=flex,
            expt_data=expt_data,
            ligand_ars=ligand_ars,
        )

        heatmap(
            symm_to_c=symm_to_c,
            results=results,
            output_dir=output_dir,
            filename=f'energy_map_flat_{flex}.pdf',
            vmin=0,
            vmax=10,
            clabel='energy (eV)',
            flex=flex,
            expt_data=expt_data,
            ligand_ars=ligand_ars,
        )

        heatmap(
            symm_to_c=symm_to_c,
            results=results,
            output_dir=output_dir,
            filename=f'shape_map_{flex}.pdf',
            vmin=0,
            vmax=2.2,
            clabel='CU-8',
            flex=flex,
            expt_data=expt_data,
            ligand_ars=ligand_ars,
        )

        scatter(
            symm_to_c=symm_to_c,
            results=results,
            output_dir=output_dir,
            filename=f'energy_{flex}.pdf',
            ylabel='energy (eV)',
            flex=flex,
        )
        scatter(
            symm_to_c=symm_to_c,
            results=results,
            output_dir=output_dir,
            filename=f'shape_{flex}.pdf',
            ylabel='CU-8',
            flex=flex,
        )

        comp_sets = {
            'ts': ('th1', 'th2', 'td', 'tl'),
            'ds': ('d31', 'd32', 'd31n', 'd32n'),
            'ss': ('s61', 's62', 's41', 's42'),
            'expt': ('d2', 'tl', 's62', 'th2', 'd32'),
        }
        for key, values in comp_sets.items():
            if flex == 'high':
                eylim = (0, 10)
            else:
                eylim = (0, 45)

            comp_scatter(
                symm_to_c=symm_to_c,
                symm_set=values,
                results=results,
                output_dir=output_dir,
                filename=f'comp_energy_{flex}_{key}.pdf',
                ylabel='energy (eV)',
                flex=flex,
                ylim=eylim,
            )
            comp_scatter(
                symm_to_c=symm_to_c,
                symm_set=values,
                results=results,
                output_dir=output_dir,
                filename=f'comp_shape_{flex}_{key}.pdf',
                ylabel='CU-8',
                flex=flex,
                ylim=(0, 2),
            )


if __name__ == "__main__":
    main()
