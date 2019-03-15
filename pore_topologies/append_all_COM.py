#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to append window and cage COM's to a CIF.

Author: Andrew Tarzia

Date Created: 19 Feb 2019
"""

import sys
from ase import Atom
from ase.io import read
from ase.visualize import view
import pywindow as pw

def convert_CIF_2_PDB(file):
    '''Convert CIF to PDB file and output for pywindow to read.

    '''
    pdb_file = file.replace('.cif', '.pdb')
    print('converting:', file, 'to', pdb_file)
    structure = read(file)
    # view(structure)
    # input()
    structure.write(pdb_file)
    print('conversion done.')
    return pdb_file, structure


def rebuild_system(file):
    '''As per example 6 in pywindow - rebuild the PDB system, output and reread.

    '''
    print('rebuilding:', file)
    molsys = pw.MolecularSystem.load_file(file)
    rebuild_molsys = molsys.rebuild_system()
    # output
    rebuild_molsys.dump_system(file.replace('.pdb', '_rebuild.pdb'),
                               include_coms=True,
                               override=True)
    print('rebuild done.')
    return rebuild_molsys


def run_analysis(rebuilt_structure, file_prefix, verbose=False):
    '''Run all desired analysis on each molecule in rebuilt structure.
        (modified version of Example6 of pywindow examples.)

    Output COM coodinates for all molecules to dictionary to add to ASE
    structure.

    '''
    result_dict = {}
    for molecule in rebuilt_structure.molecules:
        print('Analysing molecule {0} out of {1}'.format(
            molecule + 1, len(rebuilt_structure.molecules)))
        mol = rebuilt_structure.molecules[molecule]
        if mol.no_of_atoms < 20:
            continue
        try:
            analysis = mol.full_analysis()
        except ValueError as e:
            print(e)
            print('------------ passed --------------')
            continue
        if verbose:
            print(analysis, '\n')
        # Each molecule can be saved separately
        mol.dump_molecule(
            file_prefix+"_{0}.pdb".format(molecule),
            include_coms=True,
            override=True)
        mol.dump_properties_json(
            file_prefix+"_{0}.json".format(molecule),
            override=True)
        # output COM, window size and COM, and COM of pore optimized
        result_dict[molecule] = (analysis['centre_of_mass'],
                                 analysis['windows'],
                                 analysis['pore_diameter_opt']['centre_of_mass'])
    # print(result_dict)
    print('analysis done.')
    return result_dict


def append_and_write(result_dict, structure, file):
    '''Append all COMs in result dict as the atoms below to the ASE structure
    and output to file.

    Window COMs: He
    cage COM: Ne
    optimized pore COM: Ar

    '''
    for molecule in result_dict:
        # check if the molecule has any windows:
        if result_dict[molecule][1]['diameters'] is None:
            continue
        # add cage COM
        cage_com = Atom(symbol='Ne',
                        position=result_dict[molecule][0])
        structure.append(cage_com)
        # add window COMs:
        for wCOM in result_dict[molecule][1]['centre_of_mass']:
            print('----------')
            print(wCOM)
            window_com = Atom(symbol='He',
                            position=wCOM)
            structure.append(window_com)
            print('----------')
        # add optimized cage pore COM
        pore_com = Atom(symbol='Ar',
                        position=result_dict[molecule][2])
        structure.append(pore_com)
    # output to CIF
    output = file.replace('.cif', '_appended.cif')
    structure.write(output, format='cif')
    output = file.replace('.cif', '_appended.pdb')
    structure.write(output)


if __name__ == "__main__":
    if (not len(sys.argv) == 2):
        print("""
Usage: append_all_COM.py CIF
    CIF: file (.cif) to analyze and add pseudo atoms to.
    """)
        sys.exit()
    else:
        file = sys.argv[1]
    pdb_file, ASE_structure = convert_CIF_2_PDB(file)
    # rebuild system
    rebuilt_structure = rebuild_system(file=pdb_file)
    # sys.exit()
    rebuilt_structure.make_modular()
    # run analysis
    COM_dict = run_analysis(rebuilt_structure,
                            file_prefix=file.replace('.cif', ''),
                            verbose=False)
    # append atoms to ASE structure as pseudo atoms and write out new CIF
    append_and_write(COM_dict, ASE_structure, file)
