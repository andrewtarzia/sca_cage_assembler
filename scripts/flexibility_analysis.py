#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to calculate flexibility measure of ligands.

Author: Andrew Tarzia

Date Created: 21 Oct 2020

"""

import numpy as np
import sys
from os.path import join, exists
from glob import glob
import matplotlib.pyplot as plt
import json

import stk
import stko

from atools import (
    build_conformers,
    calculate_molecule_planarity,
    colors_i_like,
    histogram_plot_N,
    scatter_plot,
    update_from_rdkit_conf,
    calculate_energy,
    read_gfnx2xtb_eyfile,
    crest_conformer_search,
)


def load_ligands(directory):

    ligands = {}
    for lig in glob(join(directory, f'quad2*_opt.mol')):
        l_name = lig.replace(directory, '').replace('_opt.mol', '')
        ligands[l_name] = stk.BuildingBlock.init_from_file(
            lig,
            functional_groups=[stk.BromoFactory()],
        )

    return ligands


def calculate_flex(molecule, name):
    """
    Calculate flexibility of molecule.

    Three approaches:
        1) Based on plane deviation of Bromine atoms in molecule.
        Assumes molecule has N functional groups with only 1 binder
        atom each.

        2) Based on properties of Crest-XTBFF generated conformer
        ensemble.

        3) Based on plane deviation of all atoms in molecule.

    """

    for fg in molecule.get_functional_groups():
        if len(list(fg.get_bonder_ids())) > 1:
            raise ValueError(
                f'{molecule} has functional groups with more'
                ' than 1 binder.'
            )

    # Crest part.
    if not exists(f'crst_{name}/crest_rotamers.xyz'):
        print(f'running XTBFF-Crest on {name}')
        new_molecule = crest_conformer_search(
            molecule=molecule,
            output_dir=f'crst_{name}',
            gfn_exec='/home/atarzia/software/xtb-6.3.1/bin/xtb',
            crest_exec='/home/atarzia/software/crest/crest',
            gfn_version=2,
            nc=3,
            opt_level='crude',
            charge=0,
            keepdir=False,
            cross=False,
            etemp=300,
            no_unpaired_e=0,
            speed_setting='squick',
            solvent=('acetonitrile', 'normal'),
        )
        new_molecule.write(f'{name}_loweconf.mol')

    # Extract some measure of conformer ensemble size.
    crest_output_file = f'{name}_flex_measure.json'
    crest_data = {}
    with open(f'crst_{name}/crest.output', 'r') as f:
        for line in f.readlines():
            # Get number of conformers.
            if ' number of unique conformers for further calc' in line:
                crest_data['no_conformers'] = (
                    int(line.rstrip().split(' ')[-1])
                )

            # Get number of rotamers.
            if 'total number unique points considered further' in line:
                crest_data['no_rotamers'] = (
                    int(line.rstrip().split(' ')[-1])
                )

    # Analyse all rotamers from CREST.
    crest_conformer_files = split_xyz_file(
        num_atoms=molecule.get_num_atoms(),
        xyz_file=f'crst_{name}/crest_conformers.xyz',
    )
    bromo_ids = [
        fg.get_bromine().get_id()
        for fg in molecule.get_functional_groups()
    ]
    crest_data['binder_plane_deviations'] = [
        calculate_molecule_planarity(
            mol=molecule.with_structure_from_file(i),
            atom_ids=bromo_ids,
        )
        for i in crest_conformer_files
    ]
    crest_data['plane_deviations'] = [
        calculate_molecule_planarity(
            mol=molecule.with_structure_from_file(i),
        )
        for i in crest_conformer_files
    ]
    print(f'{name} has {len(crest_conformer_files)} conformers')

    pair1_cents = [
        molecule.with_structure_from_file(i).get_centroid(
            atom_ids=tuple(la_pairs[0])
        )
        for i in crest_conformer_files
    ]
    pair2_cents = [
        molecule.with_structure_from_file(i).get_centroid(
            atom_ids=tuple(la_pairs[1])
        )
        for i in crest_conformer_files
    ]

    la_dist = [
        np.linalg.norm(i-j)
        for i, j in zip(pair1_cents, pair2_cents)
    ]
    crest_data['lap_dist'] = la_dist

    with open(crest_output_file, 'w') as f:
        json.dump(crest_data, f)

    # Plane dev part.
    bpd_json_file = f'{name}_planedev_dist.json'
    pd_json_file = f'{name}_AAplanedev_dist.json'
    ey_json_file = f'{name}_energy_dist.json'
    ela_json_file = f'{name}_ela_dist.json'
    if (
        exists(pd_json_file) and
        exists(bpd_json_file) and
        exists(ey_json_file) and
        exists(ela_json_file)
    ):
        with open(bpd_json_file, 'r') as f:
            binder_plane_deviations = json.load(f)
        with open(pd_json_file, 'r') as f:
            plane_deviations = json.load(f)
        with open(ey_json_file, 'r') as f:
            xtb_energies = json.load(f)
        with open(ela_json_file, 'r') as f:
            ela_dist = json.load(f)
    else:

        plane_deviations = []
        binder_plane_deviations = []
        xtb_energies = []
        ela_dist = []
        cids, confs = build_conformers(
            mol=molecule,
            N=200,
            ETKDG_version='v3'
        )

        new_molecule = molecule.clone()
        for cid in cids:
            # Update stk_mol to conformer geometry.
            new_molecule = update_from_rdkit_conf(
                stk_mol=new_molecule,
                rdk_mol=confs,
                conf_id=cid
            )

            binder_plane_deviations.append(
                calculate_molecule_planarity(
                    mol=new_molecule,
                    atom_ids=bromo_ids,
                )
            )
            new_molecule.write(f'c_{name}_{cid}_ey.mol')
            plane_deviations.append(
                calculate_molecule_planarity(new_molecule)
            )
            # Calculate conformer energy.
            calculate_energy(
                name=f'{name}_{cid}',
                mol=new_molecule,
                gfn_exec=(
                    '/home/atarzia/software/xtb-6.3.1/bin/xtb'
                ),
                ey_file=f'{name}_{cid}_ey.ey',
                charge=0,
                no_unpaired_e=0,
                solvent=None

            pair1_cents = new_molecule.get_centroid(
                atom_ids=tuple(la_pairs[0])
            )
            xtb_energies.append(
                read_gfnx2xtb_eyfile(f'{name}_{cid}_ey.ey')
            pair2_cents = new_molecule.get_centroid(
                atom_ids=tuple(la_pairs[1])
            )
            ela_dist.append(np.linalg.norm(pair1_cents-pair2_cents))

        with open(bpd_json_file, 'w') as f:
            json.dump(binder_plane_deviations, f)
        with open(pd_json_file, 'w') as f:
            json.dump(plane_deviations, f)
        with open(ey_json_file, 'w') as f:
            json.dump(xtb_energies, f)
        with open(ela_json_file, 'w') as f:
            json.dump(ela_dist, f)

    return (
        binder_plane_deviations,
        plane_deviations,
        crest_data,
        xtb_energies,
        la_dist,
        ela_dist,
    )


def plot_lap(measures, name, crest=False):
    # Can assume the first one in the list of measures is the lowest
    # energy conformer.
    fig, ax = histogram_plot_N(
        Y=[i-measures[0] for i in measures],
        X_range=(-1.5, 1.5),
        width=0.05,
        alpha=1.0,
        color=colors_i_like()[1],
        edgecolor=colors_i_like()[1],
        xtitle=r'long axis deviation [$\mathrm{\AA}$]',
        N=1
    )
    range = abs(max(measures) - min(measures))
    ax.set_title(f'range = {round(range, 2)}', fontsize=16)
    fig.tight_layout()
    if crest:
        filename = f'{name}_lapC_dist.pdf'
    else:
        filename = f'{name}_lap_dist.pdf'
    fig.savefig(
        filename,
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_pd(measures, name, crest=False):

    fig, ax = histogram_plot_N(
        Y=measures,
        X_range=(0, 200),
        width=4,
        alpha=1.0,
        color=colors_i_like()[1],
        edgecolor=colors_i_like()[1],
        xtitle=r'plane deviation [$\mathrm{\AA}$]',
        N=1
    )
    range = abs(max(measures) - min(measures))
    ax.set_title(f'range = {round(range, 2)}', fontsize=16)
    fig.tight_layout()
    if crest:
        filename = f'{name}_AAplanedevC_dist.pdf'
    else:
        filename = f'{name}_AAplanedev_dist.pdf'
    fig.savefig(
        filename,
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_bpd(measures, name, crest=False):

    fig, ax = histogram_plot_N(
        Y=measures,
        X_range=(0, 30),
        width=0.2,
        alpha=1.0,
        color=colors_i_like()[1],
        edgecolor=colors_i_like()[1],
        xtitle=r'all atom plane deviation [$\mathrm{\AA}$]',
        N=1
    )
    range = abs(max(measures) - min(measures))
    ax.set_title(f'range = {round(range, 2)}', fontsize=16)
    fig.tight_layout()
    if crest:
        filename = f'{name}_planedevC_dist.pdf'
    else:
        filename = f'{name}_planedev_dist.pdf'
    fig.savefig(
        filename,
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_pd_v_ey(energies, measures, name):

    fig, ax = scatter_plot(
        X=measures,
        Y=[i-min(energies) for i in energies],
        alpha=1.0,
        xlim=None,
        ylim=None,
        c=colors_i_like()[1],
        s=60,
        edgecolors='k',
        xtitle=r'all atom plane deviation [$\mathrm{\AA}$]',
        ytitle='rel. xtb energy [kJ/mol]',
    )
    fig.tight_layout()
    fig.savefig(f'{name}_AAPDvE.pdf', dpi=720, bbox_inches='tight')
    plt.close()


def plot_bpd_v_ey(energies, measures, name):

    fig, ax = scatter_plot(
        X=measures,
        Y=[i-min(energies) for i in energies],
        alpha=1.0,
        xlim=None,
        ylim=None,
        c=colors_i_like()[1],
        s=60,
        edgecolors='k',
        xtitle=r'binder atom plane deviation [$\mathrm{\AA}$]',
        ytitle='rel. xtb energy [kJ/mol]',
    )
    fig.tight_layout()
    fig.savefig(f'{name}_plPDvE.pdf', dpi=720, bbox_inches='tight')
    plt.close()


def plot_pd_vs_crest(pd, crest):

    fig, ax = scatter_plot(
        X=pd,
        Y=crest,
        xtitle=r'std.dev. all atom PD [$\mathrm{\AA}$]',
        ytitle=r'num. crest rotamers',
        xlim=(0, None),
        ylim=(0, None),
        c=colors_i_like()[1],
        edgecolors='k',
    )

    fig.tight_layout()
    fig.savefig(
        f'pd_vs_rotamers.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_bpd_vs_crest(bpd, crest):

    fig, ax = scatter_plot(
        X=bpd,
        Y=crest,
        xtitle=r'std.dev. binder atom PD [$\mathrm{\AA}$]',
        ytitle=r'num. crest rotamers',
        xlim=(0, None),
        ylim=(0, None),
        c=colors_i_like()[1],
        edgecolors='k',
    )

    fig.tight_layout()
    fig.savefig(
        f'bpd_vs_rotamers.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def get_long_axis_atoms(ligands):

    long_axis_atom_pairs = {}
    for ligand in ligands:
        molecule = ligands[ligand]
        binder_atom_ids = [
            list(fg.get_bonder_ids())
            for fg in molecule.get_functional_groups()
        ]
        binder_atom_dists = sorted(
            [
                (idx1, idx2, get_atom_distance(
                    molecule,
                    idx1,
                    idx2
                ))
                for idx1, idx2 in combinations(binder_atom_ids, r=2)
            ],
            key=lambda a: a[2]
        )
        # Can assume the ordering of the binder atom distances:
        # 0, 1: short vectors
        # 2, 3: long vectors
        # 4, 5: diagonal vectors
        # This fails when the molecule is not sufficiently anisotropic,
        # at which point it will not matter.
        short_vector_fg_1 = (
            binder_atom_dists[0][0], binder_atom_dists[0][1]
        )
        short_vector_fg_2 = (
            (binder_atom_dists[1][0], binder_atom_dists[1][1])
            if (
                binder_atom_dists[1][0] not in short_vector_fg_1 and
                binder_atom_dists[1][1] not in short_vector_fg_1
            ) else
            (binder_atom_dists[2][0], binder_atom_dists[2][1])
        )

        pairs = (
            [i[0] for i in short_vector_fg_1],
            [i[0] for i in short_vector_fg_2]
        )
        long_axis_atom_pairs[ligand] = pairs

    return long_axis_atom_pairs


def main():
    first_line = (
        'Usage: flexibility_analysis.py '
        'lig_directory'
    )
    if (not len(sys.argv) == 2):
        print(f"""
{first_line}

    ligand_directory : (str)
        Directory with required ligand structures.

    """)
        sys.exit()
    else:
        ligand_directory = sys.argv[1]

    # Load in each ligand structure.
    ligands = load_ligands(ligand_directory)
    long_axis_atom_pairs = get_long_axis_atoms(ligands)
    print(long_axis_atom_pairs)

    crest_simplified = []
    bpd_simplified = []
    pd_simplified = []
    for lig in sorted(ligands):
        lig_structure = ligands[lig]
        bplane_devs, plane_devs, crest_data, xtb_eys, la_dist, ela_dist = (
            calculate_flex(
                molecule=lig_structure,
                name=lig,
                la_pairs=long_axis_atom_pairs[lig],
            )
        )

        print(
            lig,
            max(bplane_devs),
            np.average(bplane_devs),
            np.std(bplane_devs),
            abs(max(bplane_devs)-min(bplane_devs)),
            sum([
                abs(i-np.mean(bplane_devs)) for i in bplane_devs
            ])/len(bplane_devs)
        )
        print('>', abs(
            max(crest_data['binder_plane_deviations'])-min(
                crest_data['binder_plane_deviations']
            ))
        )
        print(
            lig,
            max(plane_devs),
            np.average(plane_devs),
            np.std(plane_devs),
            abs(max(plane_devs)-min(plane_devs)),
            sum([
                abs(i-np.mean(plane_devs)) for i in plane_devs
            ])/len(plane_devs)
        )
        print('>', abs(
            max(crest_data['plane_deviations'])-min(
                crest_data['plane_deviations']
            ))
        )
        print(
            ':::',
            abs(
                max([i-la_dist[0] for i in la_dist])
                - min([i-la_dist[0] for i in la_dist])
            )
        )
        print(
            ':::',
            abs(
                max([i-ela_dist[0] for i in ela_dist])
                - min([i-ela_dist[0] for i in ela_dist])
            )
        )
        plot_lap(la_dist, lig, crest=True)
        plot_bpd(bplane_devs, lig)
        plot_pd(plane_devs, lig)
        plot_bpd(
            crest_data['binder_plane_deviations'],
            lig,
            crest=True
        )
        plot_pd(crest_data['plane_deviations'], lig, crest=True)
        plot_bpd_v_ey(xtb_eys, bplane_devs, lig)
        plot_pd_v_ey(xtb_eys, plane_devs, lig)
        print(
            lig,
            crest_data['no_rotamers'],
            crest_data['no_conformers']
        )
        print('---')
        pd_simplified.append(abs(max(plane_devs)-min(plane_devs)))
        bpd_simplified.append(abs(max(bplane_devs)-min(bplane_devs)))
        crest_simplified.append(crest_data['no_rotamers'])

    plot_bpd_vs_crest(crest=crest_simplified, bpd=bpd_simplified)
    plot_pd_vs_crest(crest=crest_simplified, pd=pd_simplified)
    sys.exit()


if __name__ == '__main__':
    main()
