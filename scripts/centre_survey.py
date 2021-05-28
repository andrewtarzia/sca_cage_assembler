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

from plotting import (
    histogram_plot_N,
    colors_i_like,
    parity_plot,
)
from utilities import calculate_sites_order_values, get_element_sites


def plot_distances(data, atomic_number):

    distance_list = list(data.DIST1)+list(data.DIST2)
    distance_list += list(data.DIST3)+list(data.DIST4)
    distance_list += list(data.DIST5)+list(data.DIST6)

    fig, ax = histogram_plot_N(
        Y=distance_list,
        X_range=(1.5, 3.0),
        width=0.025,
        alpha=1.0,
        color=colors_i_like()[2],
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
    angle_list += list(data.ANG5)+list(data.ANG6)
    angle_list += list(data.ANG7)+list(data.ANG8)
    angle_list += list(data.ANG9)+list(data.ANG10)
    angle_list += list(data.ANG11)+list(data.ANG12)

    angle_list1 = [i for i in angle_list if i < 120]
    angle_list2 = [i for i in angle_list if i > 120]

    fig, ax = histogram_plot_N(
        Y=[angle_list1, angle_list2],
        X_range=(0, 200),
        width=5,
        alpha=[1.0, 1.0],
        color=[colors_i_like()[0], colors_i_like()[1]],
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
    order_values = calculate_sites_order_values(
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


def plot_orderparams(struct_dir, atomic_number, CN):

    # Just in case.
    if struct_dir[-1] == '/':
        file_list = sorted(glob.glob(f'{struct_dir}*.cif'))
    else:
        file_list = sorted(glob.glob(f'{struct_dir}/*.cif'))
    # If calculations have been done, then do not repeat.
    op_res_file = f'op_results{atomic_number}.json'
    if exists(op_res_file):
        with open(op_res_file, 'r') as f:
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
        exceptions = [
            'ACOMAG_extracted.cif'
        ]
        for file in file_list:
            print(file)
            if struct_dir[-1] == '/':
                clean_file = file.replace(struct_dir, '')
            else:
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
            metal_site_ids = get_element_sites(
                pmg_struct,
                atomic_no=atomic_number
            )
            # Define neighbouring N atoms.
            v = CrystalNN()
            centre_mols = {}
            for msite in metal_site_ids:
                crossing = False
                non_nitrogrens = False
                nninfo = v.get_nn_info(pmg_struct, n=msite)
                Nsites = []
                for i in nninfo:
                    # Remove metal sites that cross boundaries.
                    if i['image'] != (0, 0, 0):
                        crossing = True
                    isite = pmg_struct.species[i['site_index']]
                    isite_symbol = isite.symbol
                    if isite_symbol != 'N':
                        non_nitrogrens = True
                        break
                    else:
                        Nsites.append(i['site_index'])
                # Skip msite if any elements != N or N_sites != CN.
                if non_nitrogrens:
                    print(
                        f'for {msite}, an element is non-Nitrogen '
                        f'{isite} - {isite_symbol}'
                    )
                    continue
                if len(Nsites) != CN:
                    print(
                        f'for {msite}, no. Nsites != {CN} - skipping'
                        f'\n{Nsites} :: '
                        f'{[pmg_struct[i] for i in Nsites]}'
                    )
                    continue
                if not crossing:
                    coords = [pmg_struct[msite].coords]
                    for site in Nsites:
                        coords.append(pmg_struct[site].coords)
                    m_symbol = pmg_struct.species[msite]
                    centre_mol = Molecule(
                        species=[m_symbol]+['N']*CN,
                        coords=coords
                    )
                    centre_mols[msite] = centre_mol

            # Calculate order parameters.
            for msite in centre_mols:
                cmol = centre_mols[msite]
                order_values = calculate_sites_order_values(
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
                results['pdsite'].append(msite)
                results['sqpl'].append(OPs['sq_plan'])
                results['oct'].append(OPs['oct'])
                results['q4'].append(OPs['q4'])
                results['q6'].append(OPs['q6'])

        with open(op_res_file, 'w') as f:
            json.dump(results, f)

    # for i, c in enumerate(results['coords']):
    #     if results['sqpl'][i] < 0.6:
    #         print(results['sqpl'][i])
    #         input(c)

    fig, ax = histogram_plot_N(
        Y=results['sqpl'],
        X_range=(0, 1),
        width=0.05,
        alpha=1.0,
        density=True,
        color=colors_i_like()[0],
        edgecolor='k',
        xtitle=r'$q_{\mathrm{sqp}}$'
    )
    fig.tight_layout()
    fig.savefig(
        f'sqpl_{atomic_number}_survey.pdf',
        dpi=720,
        bbox_inches='tight'
    )

    fig, ax = histogram_plot_N(
        Y=results['oct'],
        X_range=(0, 1),
        width=0.05,
        alpha=1.0,
        density=True,
        color=colors_i_like()[1],
        edgecolor='k',
        xtitle=r'$q_{\mathrm{oct}}$'
    )
    fig.tight_layout()
    fig.savefig(
        f'oct_{atomic_number}_survey.pdf',
        dpi=720,
        bbox_inches='tight'
    )

    fig, ax = histogram_plot_N(
        Y=results['q4'],
        X_range=(0, 1),
        width=0.05,
        alpha=1.0,
        density=True,
        color=colors_i_like()[5],
        edgecolor='k',
        xtitle=r'$q_4$'
    )
    fig.tight_layout()
    fig.savefig(
        f'q4_{atomic_number}_survey.pdf',
        dpi=720,
        bbox_inches='tight'
    )

    fig, ax = histogram_plot_N(
        Y=results['q6'],
        X_range=(0, 1),
        width=0.05,
        alpha=1.0,
        density=True,
        color=colors_i_like()[2],
        edgecolor='k',
        xtitle=r'$q_6$'
    )
    fig.tight_layout()
    fig.savefig(
        f'q6_{atomic_number}_survey.pdf',
        dpi=720,
        bbox_inches='tight'
    )

    fig, ax = parity_plot(
        X=results['sqpl'],
        Y=results['oct'],
        lim=(0, 1),
        c=colors_i_like()[2],
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
    if (not len(sys.argv) == 5):
        print("""
    Usage: centre_survey.py file struct_dir atomic_number
        file (str) - csv file to analyze
        struct_dir (str) - dir with struct files to analyze
        atomic_number (str) - atomic number of metal to analyse
        CN (str) - integer number of nitrogen neighbours of metal
        """)
        sys.exit()
    else:
        file = sys.argv[1]
        struct_dir = sys.argv[2]
        atomic_number = int(sys.argv[3])
        CN = int(sys.argv[4])

    data = pd.read_csv(file)
    print(data.columns)

    plot_orderparams(struct_dir, atomic_number, CN)
    plot_distances(data, atomic_number)
    plot_angles(data, atomic_number)


if __name__ == "__main__":
    main()
