#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions that are useful for pywindow usage

Author: Andrew Tarzia

Date Created: 15 Mar 2019
"""

from rdkit.Chem import AllChem as Chem
from ase.atoms import Atom
import pywindow as pw


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


def analyze_cage(cage, propfile=None, structfile=None, include_coms=True):
    '''Analyze cage already loaded into pywindow.

    '''
    # Perform full pyWindow analysis
    cage.full_analysis()
    # Dump pyWindow properties into JSON and cage into xyz
    if propfile is not None:
        cage.dump_properties_json(propfile)
    if structfile is not None:
        cage.dump_molecule(structfile, include_coms=include_coms)


def analyze_cage_from_MOL(file, prop_file, mole_file, include_coms=True):
    '''Run all desired analysis on a single built cage molecule.

    Output cage with COM atoms and properties to JSON.

    '''
    # Import optimised cage into pyWindow, via RDkit mol file
    cage_rd = Chem.MolFromMolFile(file)
    cage_sys = pw.MolecularSystem.load_rdkit_mol(cage_rd)
    cage_mol = cage_sys.system_to_molecule()
    analyze_cage(cage=cage_mol, propfile=prop_file,
                 structfile=mole_file, include_coms=include_coms)


def analyze_rebuilt(rebuilt_structure, file_prefix, atom_limit,
                    include_coms=True, verbose=False):
    '''Run all desired analysis on each molecule in rebuilt structure.
        (modified version of Example6 of pywindow examples.)

    Keyword Arguments:
        rebuilt_structure (Molecule) - Pywindow rebuilt unitcell with cage
            molecules
        file_prefix (str) - file naming convention
        atom_limit (int) - number of atoms used as cutoff for analysis
        include_coms (bool) - whether output PDB includes window COMs
        verbose (bool) - [Default = False]

    Returns:
        result_dict (dictionary) - dictionary of window information for all
            cages

    '''
    result_dict = {}
    for molecule in rebuilt_structure.molecules:
        print('Analysing molecule {0} out of {1}'.format(
            molecule + 1, len(rebuilt_structure.molecules)))
        mol = rebuilt_structure.molecules[molecule]
        if mol.no_of_atoms < atom_limit:
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
            file_prefix + "_{0}.pdb".format(molecule),
            include_coms=include_coms,
            override=True)
        mol.dump_properties_json(
            file_prefix + "_{0}.json".format(molecule),
            override=True)
        # output COM, window size and COM, and COM of pore optimized
        result_dict[molecule] = (analysis['centre_of_mass'],
                                 analysis['windows'],
                                 analysis['pore_diameter_opt']['centre_of_mass'])
    # print(result_dict)
    print('analysis done.')
    return result_dict


def append_and_write_COMs(result_dict, structure, file):
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
