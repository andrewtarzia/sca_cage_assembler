#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to analyze the formation properties of a list of built cages based on
energies and properties calculated into CSV files.

Author: Andrew Tarzia

Date Created: 09 Apr 2019

"""

from pandas import read_csv
import sys
sys.path.insert(0, '/home/atarzia/thesource/')
from calculations import calc_formation_energy
from stk_functions import topo_2_stoich


def main():
    if (not len(sys.argv) == 6):
        print("""
Usage: analyze_molecule_list.py cage_csv alde_csv amine_csv other_csv energy
    cage_csv (str) - CSV file with cage properties
    alde_csv (str) - CSV file with aldehyde precursor properties
    amine_csv (str) - CSV file with amine precursor properties
    other_csv (str) - CSV file with other molecule (water) properties
    energy (str) - FE for free energy, TE for total energy
    """)
        sys.exit()
    else:
        cage_csv = sys.argv[1]
        alde_csv = sys.argv[2]
        amin_csv = sys.argv[3]
        othe_csv = sys.argv[4]
        energy = sys.argv[5]

    if energy == 'FE':
        used_energy = 'FE_au'
    elif energy == 'TE':
        used_energy = 'TE_au'

    print('-----------------------------------------')
    print('doing all sorts of analysis using {}.'.format(used_energy))
    print('-----------------------------------------')

    cage_db = read_csv(cage_csv)
    alde_db = read_csv(alde_csv)
    amin_db = read_csv(amin_csv)
    othe_db = read_csv(othe_csv)

    # properties
    prop = {}

    # iterate over cages
    total_cages = len(cage_db)
    total_cages_calc_done = len(cage_db[cage_db.calc_complete == 0])
    percent_calc_done = total_cages_calc_done / total_cages
    print('{} of {} ({:2.2%}) cages have complete calculations'.format(total_cages_calc_done, total_cages, percent_calc_done))
    for i, row in cage_db.iterrows():
        # focus on one aldehyde for now
        if row.aldehyde != 1:
            continue
        # is calculation done?
        if row.calc_complete == 1:
            print('{} calculation incomplete.'.format(row.file))
            continue
        # calculate cage formation energy
        # collect energies of all components in their stoichiometries
        prod_energies = [row[used_energy]]
        # add N waters, N == no imines
        print(prod_energies, row.no_imines)
        water_energy = [float(othe_db[othe_db.file == 'water.xyz'][used_energy])]
        print(water_energy)
        prod_energies += water_energy * row.no_imines
        print(prod_energies)
        react_energies = []
        print(react_energies, row.aldehyde, row.amine, row.topology)
        # topo_2_stoich outputs the stoich in order of decreasing connectivity
        print('assuming aldehyde is 2 connected, amine is 3 connected!!')
        amin_stoich, alde_stoich = topo_2_stoich(row.topology)
        print(alde_stoich, amin_stoich)
        alde_ey = [float(alde_db[alde_db.file == str(row.aldehyde) + '.xyz'][used_energy])]
        print(alde_ey)
        react_energies += alde_ey * alde_stoich
        print(react_energies)
        amin_ey = [float(amin_db[amin_db.file == str(row.amine) + '.xyz'][used_energy])]
        print(amin_ey)
        react_energies += amin_ey * amin_stoich
        print(react_energies)
        c_fe = calc_formation_energy(prod_energies, react_energies)
        c_fe_perimine = c_fe / row.no_imines
        print(c_fe)
        print(c_fe_perimine)
        prop[row.file] = (c_fe, c_fe_perimine, row.shape_persist)
        input()

    print(prop)
    for i in prop:
        if '_B_' in i:
            print(i, prop[i][1] * 2625.50, prop[i][2])
    sys.exit()


if __name__ == "__main__":
    main()
