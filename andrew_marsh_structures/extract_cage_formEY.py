#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to calculate the formation energies of a list of cages.

Author: Andrew Tarzia

Date Created: 03 May 2019
"""

import sys
import os
import glob
import pandas as pd
sys.path.insert(0, '/home/atarzia/thesource/')
import stk_f
import plotting


def get_energy_from_DF(DF, name):
    '''Return the float of the FE ('FE (au)')) from a pandas DataFrame.

    '''
    rowDF = DF[DF['NAMES'] == name].iloc[0]
    value = rowDF['FE (au)']
    return float(value)


def main():
    """Run script.

    """
    if (not len(sys.argv) == 4):
        print("""
Usage: get_all_struct_energies.py output_file suffix
    output_file (str): file to output results
    cage_energy (str): file with cage energies
    prec_energy (str): file with energies of all required precursors
""")
        sys.exit()
    else:
        output_file = sys.argv[1]
        cage_file = sys.argv[2]
        prec_file = sys.argv[3]

    if os.path.isfile(output_file):
        # read in data
        DATA = pd.read_csv(output_file)
    else:
        DATA = pd.DataFrame(columns=['file', 'cagename', 'bb1', 'bb2', 'topo', 'FE (au)'])

    # read in energy files
    cage_energies = pd.read_csv(cage_file)
    prec_energies = pd.read_csv(prec_file)
    print(cage_energies)
    print('---')
    print(prec_energies)

    done_files = list(DATA['file'])
    cage_opt_files = glob.glob('*_opt.mol')

    # remove the cages for which the formation energies have been output
    cage_opt_files = [i for i in cage_opt_files if i not in done_files]
    # get the formation energies of the remaining, while updating the output file
    for file in cage_opt_files:
        print(file)
        # naming
        name = file.replace('_opt.mol', '')
        bb1 = name.split('_')[0]  # alde
        bb2 = '_'.join(name.split('_')[1:3])  # amine
        topo = name.split('_')[3]
        # get cage energy
        cage_energy = get_energy_from_DF(DF=prec_energies, name=name)
        print('c', cage_energy)
        # get water energy * noimines formed
        noimines = stk_f.topo_2_property(topology=topo,
                                                 property='noimines')
        water_energy = get_energy_from_DF(DF=prec_energies, name='water') * noimines
        print('w', water_energy, water_energy / noimines)
        # get bb1 energy * stoich
        bb1_stoich, bb2_stoich = stk_f.topo_2_property(topology=topo,
                                                               property='stoich')
        bb1_energy = get_energy_from_DF(DF=prec_energies, name=bb1) * bb1_stoich
        print('bb1', bb1_energy, bb1_energy / bb1_stoich, bb1_stoich)
        # get bb2 energy * stoich
        bb2_energy = get_energy_from_DF(DF=prec_energies, name=bb2) * bb2_stoich
        print('bb2', bb2_energy, bb2_energy / bb2_stoich, bb2_stoich)
        # calc formation energy
        form_energy = (cage_energy + water_energy) - (bb1_energy + bb2_energy)
        print(form_energy)
        DATA = DATA.append({'file': file,
                            'cagename': name,
                            'bb1': bb1,
                            'bb2': bb2,
                            'topo': topo,
                            'FE (au)': form_energy},
                           ignore_index=True)
        print(DATA)
        sys.exit()
        # write dataFrame
        DATA.to_csv(output_file, index=False)

    # plot histogram of formation energies
    fig, ax = plotting.histogram_plot_1(Y=DATA['FE (au)'],
                                        X_range=(-1000, 1000),
                                        width=10,
                                        alpha=1,
                                        color='mediumvioletred',
                                        edgecolor='k',
                                        xtitle='formation energy [a.u.]')
    fig.tight_layout()
    fig.savefig('cage_formEY_GFN.pdf', dpi=720,
                bbox_inches='tight')


if __name__ == "__main__":
    main()
