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

    with open(crest_output_file, 'w') as f:
        json.dump(crest_data, f)

    # Plane dev part.
    bpd_json_file = f'{name}_planedev_dist.json'
    pd_json_file = f'{name}_AAplanedev_dist.json'
    ey_json_file = f'{name}_energy_dist.json'
    if (
        exists(pd_json_file) and
        exists(bpd_json_file) and
        exists(ey_json_file)
    ):
        with open(bpd_json_file, 'r') as f:
            binder_plane_deviations = json.load(f)
        with open(pd_json_file, 'r') as f:
            plane_deviations = json.load(f)
        with open(ey_json_file, 'r') as f:
            xtb_energies = json.load(f)
    else:
        bromo_ids = [
            fg.get_bromine().get_id()
            for fg in molecule.get_functional_groups()
        ]

        plane_deviations = []
        binder_plane_deviations = []
        xtb_energies = []
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
            )
            xtb_energies.append(
                read_gfnx2xtb_eyfile(f'{name}_{cid}_ey.ey')
            )
        with open(bpd_json_file, 'w') as f:
            json.dump(binder_plane_deviations, f)
        with open(pd_json_file, 'w') as f:
            json.dump(plane_deviations, f)
        with open(ey_json_file, 'w') as f:
            json.dump(xtb_energies, f)

    return (
        binder_plane_deviations,
        plane_deviations,
        crest_data,
        xtb_energies,
    )


def plot_pd(measures, name):

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
    fig.tight_layout()
    fig.savefig(
        f'{name}_AAplanedev_dist.pdf',
        dpi=720,
        bbox_inches='tight'
    )


def plot_bpd(measures, name):

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
    fig.tight_layout()
    fig.savefig(
        f'{name}_planedev_dist.pdf',
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
        xlim=(0, 50),
        ylim=(0, 4000),
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
        xlim=(0, 30),
        ylim=(0, 4000),
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
    crest_simplified = []
    bpd_simplified = []
    pd_simplified = []
    for lig in sorted(ligands):
        lig_structure = ligands[lig]
        bplane_devs, plane_devs, crest_data, xtb_eys = (
            calculate_flex(
                molecule=lig_structure,
                name=lig,
            )
        )
        print(
            lig,
            max(bplane_devs),
            np.average(bplane_devs),
            np.std(bplane_devs)
        )
        print(
            lig,
            max(plane_devs),
            np.average(plane_devs),
            np.std(plane_devs)
        )
        plot_bpd(bplane_devs, lig)
        plot_pd(plane_devs, lig)
        plot_bpd_v_ey(xtb_eys, bplane_devs, lig)
        plot_pd_v_ey(xtb_eys, plane_devs, lig)
        print(lig, crest_data)
        print('---')
        pd_simplified.append(np.std(plane_devs))
        bpd_simplified.append(np.std(bplane_devs))
        crest_simplified.append(crest_data['no_rotamers'])

    plot_bpd_vs_crest(crest=crest_simplified, bpd=bpd_simplified)
    plot_pd_vs_crest(crest=crest_simplified, pd=pd_simplified)
    sys.exit()


if __name__ == '__main__':
    main()
