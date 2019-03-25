#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to get RMSD between two structures.

Author: Andrew Tarzia

Date Created: 08 Mar 2019

"""

import sys
import os
import glob
import json
import pywindow as pw
sys.path.insert(0, '/home/atarzia/thesource/')
import conversion
import plotting
import pywindow_functions
import stk_functions


def main():
    if (not len(sys.argv) == 1):
        print("""
Usage: test_shape.py
    """)
        sys.exit()
    else:
        pass

    # get CIFs of interest -- cleaned if it exists, otherwise original
    list_of_cifs = [i for i in glob.glob('*.cif') if '_cleaned' not in i]
    list_of_cleaned_cifs = glob.glob('*_cleaned.cif')
    calculation_list = [i for i in list_of_cifs
                        if i.replace('.cif', '_cleaned.cif') not in list_of_cleaned_cifs]
    for i in list_of_cleaned_cifs:
        calculation_list.append(i)
    for calc in calculation_list:
        print(calc)
        pre_op = calc.replace('.cif', '_preop')
        pdb_file = calc.replace('.cif', '.pdb')
        print(pdb_file, pre_op)
        if os.path.isfile(pdb_file) is False:
            pdb_file, _ = conversion.convert_CIF_2_PDB(calc)
            del _  # we don't need the ASE structure in this case
        # rebuild system
        rebuilt_structure = pywindow_functions.rebuild_system(file=pdb_file)
        rebuilt_structure.make_modular()
        # run analysis on rebuilt system (extracts all cages)
        _ = pywindow_functions.analyze_rebuilt(rebuilt_structure,
                                               file_prefix=pre_op,
                                               atom_limit=20,
                                               include_coms=False,
                                               verbose=False)
        del _  # not needed again
        # determine independant cages based on pore diameters
        indep_cages = {}
        for js in glob.glob(pre_op+'*.json'):
            with open(js, 'r') as f:
                data = json.load(f)
            if data['pore_diameter_opt']['diameter'] > 0:
                ID = js.rstrip('.json')
                indep_cages[ID] = data['pore_diameter_opt']['diameter']
        print(indep_cages)
        # check shape persistency of each independant cage
        cage_output = {}
        for ID in indep_cages:
            newID = ID.replace('preop', 'postop')
            # run MD with OPLS if output file does not exist
            if os.path.isfile(newID+'.pdb') is False:
                print('doing optimization of cage:', ID)
                stk_functions.optimize_structunit(
                        infile=ID+'.pdb',
                        outfile=newID+'.pdb',
                        exec='/home/atarzia/software/schrodinger_install',
                        method='OPLS',
                        settings=None,
                        md=None)
                print('done')
            # analyze optimized cage with pyWindow and output to JSON
            if os.path.isfile(newID+'.json') is False:
                cagesys = pw.MolecularSystem.load_file(newID+'.pdb')
                cage = cagesys.system_to_molecule()
                pywindow_functions.analyze_cage(cage=cage,
                                                propfile=newID+'.json',
                                                structfile=None)
            with open(newID+'.json', 'r') as f:
                data = json.load(f)
            new_pore_diam = data['pore_diameter_opt']['diameter']
            cage_output[newID] = (indep_cages[ID], new_pore_diam)
            print('-----------------------------------------')
            print(ID, newID)
            print('preop diameter:', indep_cages[ID])
            print('postop diameter:', new_pore_diam)
            print('-----------------------------------------')
        # output plot for each CIF
        plotx = [cage_output[i][0] for i in cage_output]
        ploty = [cage_output[i][1] for i in cage_output]
        plotting.parity_plot(X=plotx, Y=ploty,
                             outfile=calc.replace('.cif', '.pdf'),
                             xtitle='XRD pore diameter [$\mathrm{\AA}$]',
                             ytitle='opt. pore diameter [$\mathrm{\AA}$]',
                             lim=(0, round(max([max(ploty), max(plotx)]))+1),
                             )


if __name__ == "__main__":
    main()
