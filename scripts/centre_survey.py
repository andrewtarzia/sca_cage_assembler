#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to analyze the properties of metal centres in the CSD.

Author: Andrew Tarzia

Date Created: 9 Mar 2020

"""

import sys
import glob
import numpy as np
import pandas as pd
from pymatgen import Molecule
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.io.cif import CifParser
from os.path import exists
import json

import atools


def plot_distances(data, atomic_number):

    distance_list = list(data.DIST1)+list(data.DIST2)
    distance_list += list(data.DIST3)+list(data.DIST4)

    fig, ax = atools.histogram_plot_N(
        Y=distance_list,
        X_range=(1.8, 2.4),
        width=0.025,
        alpha=1.0,
        color=atools.colors_i_like()[2],
        edgecolor='k',
        xtitle=r'N-M bond distance [$\mathrm{\AA}$]'
    )
    fig.tight_layout()
    fig.savefig(
        f'N_M{atomic_number}_bond_survey.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    print('bond lengths:')
    print(f'Mean = {np.mean(distance_list)} Angstrom')
    print(f'Standard Deviation = {np.std(distance_list)} Angstrom')
    print('-------------------------------------------------------')


def plot_angles(data, atomic_number):

    angle_list = list(data.ANG1)+list(data.ANG2)
    angle_list += list(data.ANG3)+list(data.ANG4)
    angle_list += list(data.V1)+list(data.V2)

    angle_list1 = [i for i in angle_list if i < 120]
    angle_list2 = [i for i in angle_list if i > 120]

    fig, ax = atools.histogram_plot_N(
        Y=[angle_list1, angle_list2],
        X_range=(0, 200),
        width=5,
        alpha=[1.0, 1.0],
        color=[atools.colors_i_like()[0], atools.colors_i_like()[1]],
        edgecolor='k',
        xtitle=r'N-M-N angles [$^{\circ}$]',
        labels=['< 120', '> 120'],
        N=2
    )
    fig.tight_layout()
    fig.savefig(
        f'N_M{atomic_number}_angle_survey.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    print('angle set 1:')
    print(f'Mean = {np.mean(angle_list1)} Degrees')
    print(f'Standard Deviation = {np.std(angle_list1)} Degrees')
    print('-------------------------------------------------------')
    print('angle set 2:')
    print(f'Mean = {np.mean(angle_list2)} Degrees')
    print(f'Standard Deviation = {np.std(angle_list2)} Degrees')
    print('-------------------------------------------------------')


def calculate_orderparams(coord_set):
    # Create pymatgen object from coord_set.
    pmg_mol = Molecule(
        species=[i[2] for i in coord_set],
        coords=[[i[3], i[4], i[5]] for i in coord_set]
    )
    sites = [0]
    neighs = [[j for j in range(1, len(coord_set))]]

    # Calculate order parameters.
    order_values = atools.calculate_sites_order_values(
        molecule=pmg_mol,
        site_idxs=sites,
        neigh_idxs=neighs
    )
    OPs = {
        'sqpl': order_values[0]['sq_plan'],
        'oct': order_values[0]['oct'],
        'q4': order_values[0]['q4'],
        'q6': order_values[0]['q6']
    }

    return OPs


def plot_orderparams(struct_dir, atomic_number):

    # Just in case.
    clean_dir = struct_dir.replace('/', '')
    file_list = sorted(glob.glob(f'{clean_dir}/*.cif'))
    # If calculations have been done, then do not repeat.
    if exists('op_results.json'):
        with open('op_results.json', 'r') as f:
            results = json.load(f)
    else:
        results = {
            'file': [],
            'pdsite': [],
            'sqpl': [],
            'oct': [],
            'q4': [],
            'q6': []
        }
        # Get neighbours and order parameters directly from CIF with
        # pymatgen.
        exceptions = {}
        for file in file_list:
            print(file)
            clean_file = file.replace(struct_dir+'/', '')
            if clean_file in exceptions:
                continue
            # Set high occupancy_tolerance to allow reading.
            pmg_struct = CifParser(
                file,
                occupancy_tolerance=100,
                site_tolerance=1e-2
            ).get_structures()[0]
            # Then merge close sites.
            print(len(pmg_struct.species))
            pmg_struct.merge_sites(tol=0.1, mode='delete')
            # Make supercell.
            pmg_struct.make_supercell([2, 2, 2])
            print(len(pmg_struct.species))
            pd_site_ids = atools.get_element_sites(
                pmg_struct,
                atomic_no=atomic_number
            )

            # Define neighbouring N atoms.
            v = CrystalNN()
            centre_mols = {}
            for pdsite in pd_site_ids:
                crossing = False
                nninfo = v.get_nn_info(pmg_struct, n=pdsite)
                Nsites = []
                for i in nninfo:
                    # Remove pd sites that cross boundaries.
                    if i['image'] != (0, 0, 0):
                        crossing = True
                    Nsites.append(i['site_index'])
                if not crossing:
                    centre_mol = Molecule(
                        species=['Pd', 'N', 'N', 'N', 'N'],
                        coords=[
                            pmg_struct[pdsite].coords,
                            pmg_struct[Nsites[0]].coords,
                            pmg_struct[Nsites[1]].coords,
                            pmg_struct[Nsites[2]].coords,
                            pmg_struct[Nsites[3]].coords
                        ]
                    )
                    centre_mols[pdsite] = centre_mol

            # Calculate order parameters.
            for pdsite in centre_mols:
                cmol = centre_mols[pdsite]
                order_values = atools.calculate_sites_order_values(
                    molecule=cmol,
                    site_idxs=[0],
                    neigh_idxs=[[1, 2, 3, 4]],
                    target_species_type=None
                )
                print(order_values)
                OPs = order_values[0]
                if OPs['sq_plan'] < 0.5:
                    print(file)
                    print(OPs)
                results['file'].append(file)
                results['pdsite'].append(pdsite)
                results['sqpl'].append(OPs['sq_plan'])
                results['oct'].append(OPs['oct'])
                results['q4'].append(OPs['q4'])
                results['q6'].append(OPs['q6'])

        with open('op_results.json', 'w') as f:
            json.dump(results, f)

    # for i, c in enumerate(results['coords']):
    #     if results['sqpl'][i] < 0.6:
    #         print(results['sqpl'][i])
    #         input(c)

    fig, ax = atools.histogram_plot_N(
        Y=results['sqpl'],
        X_range=(0, 1),
        width=0.05,
        alpha=1.0,
        density=True,
        color=atools.colors_i_like()[0],
        edgecolor='k',
        xtitle=r'$q_{\mathrm{sqp}}$'
    )
    fig.tight_layout()
    fig.savefig(
        f'sqpl_{atomic_number}_survey.pdf',
        dpi=720,
        bbox_inches='tight'
    )

    fig, ax = atools.histogram_plot_N(
        Y=results['oct'],
        X_range=(0, 1),
        width=0.05,
        alpha=1.0,
        density=True,
        color=atools.colors_i_like()[1],
        edgecolor='k',
        xtitle=r'$q_{\mathrm{oct}}$'
    )
    fig.tight_layout()
    fig.savefig(
        f'oct_{atomic_number}_survey.pdf',
        dpi=720,
        bbox_inches='tight'
    )

    fig, ax = atools.histogram_plot_N(
        Y=results['q4'],
        X_range=(0, 1),
        width=0.05,
        alpha=1.0,
        density=True,
        color=atools.colors_i_like()[5],
        edgecolor='k',
        xtitle=r'$q_4$'
    )
    fig.tight_layout()
    fig.savefig(
        f'q4_{atomic_number}_survey.pdf',
        dpi=720,
        bbox_inches='tight'
    )

    fig, ax = atools.histogram_plot_N(
        Y=results['q6'],
        X_range=(0, 1),
        width=0.05,
        alpha=1.0,
        density=True,
        color=atools.colors_i_like()[2],
        edgecolor='k',
        xtitle=r'$q_6$'
    )
    fig.tight_layout()
    fig.savefig(
        f'q6_{atomic_number}_survey.pdf',
        dpi=720,
        bbox_inches='tight'
    )

    fig, ax = atools.parity_plot(
        X=results['sqpl'],
        Y=results['oct'],
        lim=(0, 1),
        c=atools.colors_i_like()[2],
        marker='o',
        xtitle=r'avg. $q_{\mathrm{oct}}$',
        ytitle=r'avg. $q_{\mathrm{sqp}}$'
    )
    fig.tight_layout()
    fig.savefig(
        f'sqpl_oct_{atomic_number}_survey.pdf',
        dpi=720,
        bbox_inches='tight'
    )


def main():
    if (not len(sys.argv) == 4):
        print("""
    Usage: centre_survey.py file struct_dir atomic_number
        file (str) - csv file to analyze
        struct_dir (str) - dir with struct files to analyze
        atomic_number (str) - atomic number of metal to analyse
        """)
        sys.exit()
    else:
        file = sys.argv[1]
        struct_dir = sys.argv[2]
        atomic_number = sys.argv[3]

    data = pd.read_csv(file)
    print(data.columns)

    plot_orderparams(struct_dir, atomic_number)
    plot_distances(data, atomic_number)
    plot_angles(data, atomic_number)


if __name__ == "__main__":
    main()
