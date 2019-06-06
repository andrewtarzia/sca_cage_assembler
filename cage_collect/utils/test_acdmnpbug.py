#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to test for the AZOBEY bug in all CIFs from POC_list_2

Author: Andrew Tarzia

Date Created: 22 May 2019

"""

import logging
import sys
import numpy as np
import ase
from collections import Counter
import scipy
import pywindow as pw
from pywindow.utilities import _FunctionError as _FunctionError
sys.path.insert(0, '/home/atarzia/thesource/')
import pywindow_f


def mod_mod(file):
    '''Run test for bug with pdb file.

    This block of code is a self-contained version of pw.utilities.discrete_molecules()

    '''
    molS = pw.MolecularSystem.load_file(file)
    system = molS.system
    supercell_333 = pw.utilities.create_supercell(
        system=system,
        supercell=[[-1, 1], [-1, 1], [-1, 1]])
    rebuild = supercell_333
    tol = 0.4

    # First we check which operation mode we use.
    #    1) Non-periodic MolecularSystem.
    #    2) Periodic MolecularSystem without rebuilding.
    #    3) Periodic Molecular system with rebuilding (supercell provided).
    if rebuild is not None:
        mode = 3
    else:
        if 'unit_cell' in system.keys():
            if system['unit_cell'].shape == (6,):
                mode = 2
            else:
                mode = 1
        elif 'lattice' in system.keys():
            if system['lattice'].shape == (3, 3):
                mode = 2
            else:
                mode = 1
        else:
            mode = 1

    # We create a list containing all atoms, theirs periodic elements and
    # coordinates. As this process is quite complicated, we need a list
    # which we will gradually be reducing.
    try:
        elements = system['elements']
        coordinates = system['coordinates']
    except KeyError:
        raise _FunctionError(
            "The 'elements' key is missing in the 'system' dictionary "
            "attribute of the MolecularSystem object. Which means, you need to"
            " decipher the forcefield based atom keys first (see manual)."
        )

    # print(f'elements: {len(elements)}')
    # print(f'coordinates: {len(coordinates)}')

    coordinates = system['coordinates']
    args = (elements, coordinates)
    adj = 0
    # If there are forcefield 'atom ids' as well we will retain them.
    if 'atom_ids' in system.keys():
        atom_ids = system['atom_ids']
        args = (elements, atom_ids, coordinates)
        adj = 1
    atom_list = pw.utilities.compose_atom_list(*args)
    atom_coor = pw.utilities.decompose_atom_list(atom_list)[1 + adj]

    print(atom_coor)

    dist_matrix = pw.utilities.euclidean_distances(atom_coor,
                                                   atom_coor)
    print(dist_matrix)
    dist_matrix = scipy.spatial.distance.pdist(atom_coor)
    print(dist_matrix)
    print(len(dist_matrix))
    print(min(dist_matrix))

    sys.exit()

    # print(f'atom_list: {len(atom_list)}')
    # print(f'atom_coor: {len(atom_coor)}')

    # Scenario 1: We load a non-periodic MolecularSystem.
    # We will not have 'unit_cell' nor 'lattice' keywords in the dictionary
    # and also we do not do any re-building.
    # Scenario 2: We load a periodic MolecularSystem. We want to only Extract
    # complete molecules that do not have been affected by the periodic
    # boundary.
    # Scenario 3: We load a periodic Molecular System. We want it to be rebuild
    # therefore, we also provide a supercell.
    # Scenarios 2 and 3 require a lattice and also their origin is at origin.
    # Scenario 1 should have the origin at the center of mass of the system.
    # EDIT 09-04-18: All origins/pseudo_origin had to be skewed towards some
    # direction (x + 0.01) so that there would be no ambiguity in periodic
    # ang highly symmetric systems where the choice of the closest atom would
    # be random from a set of equally far choices - bug found in the testing
    # this way rebuild system should always look the same from the same input
    # and on different machines.
    if mode == 2 or mode == 3:
        # print(f'Scen 2,3')
        # Scenarios 2 or 3.
        origin = np.array([0.01, 0., 0.])
        if 'lattice' not in system.keys():
            matrix = pw.utilities.unit_cell_to_lattice_array(system['unit_cell'])
        else:
            matrix = system['lattice']
        # print(f'matrix: {matrix}')
        pseudo_origin_frac = np.array([0.26, 0.25, 0.25])
        pseudo_origin = pw.utilities.cartisian_from_fractional(pseudo_origin_frac, matrix)
        # If a supercell is also provided that encloses the unit cell for the
        # reconstruction of the molecules through the periodic boundary.
        if rebuild is not None:
            selements = rebuild['elements']
            sids = rebuild['atom_ids']
            scoordinates = rebuild['coordinates']
            satom_list = pw.utilities.compose_atom_list(selements, sids, scoordinates)
            satom_coor = pw.utilities.decompose_atom_list(satom_list)[1 + adj]
            # print(f'satom_list: {len(satom_list)}')
            # print(f'satom_coor: {len(satom_coor)}')
        # There is one more step. We need to sort out for all the
        # reconstructed molecules, which are the ones that belong to the
        # unit cell. As we did the reconstruction to every chunk in the unit
        # cell we have now some molecules that belong to neighbouring cells.
        # The screening is simple. If the COM of a molecule translated to
        # fractional coordinates (so that it works for parallelpiped) is
        # within the unit cell boundaries <0, 1> then it's it. There is
        # an exception, for the trajectories, very often the unit cell
        # is centered at origin. Therefore we need to use <-0.5, 0.5>
        # boundary. We will simply decide which is the case by calculating
        # the centre of mass of the whole system.
        system_com = pw.utilities.center_of_mass(elements, coordinates)
        if np.allclose(system_com, origin, atol=1e-00):
            boundary = np.array([-0.5, 0.5])
        else:
            boundary = np.array([0., 1.])
        # print(f'boundary: {boundary}')
    else:
        # Scenario 1.
        pseudo_origin = pw.utilities.center_of_mass(
            elements, coordinates) + np.array([0.01, 0., 0.])
    # Here the final discrete molecules will be stored.
    molecules = []
    bools = []
    # print(f'no molecules: {len(molecules)}')
    # Exceptions. Usually end-point atoms that create single bonds or
    # just a separate atoms in the system.
    exceptions = ['H', 'CL', 'BR', 'F', 'HE', 'AR', 'NE', 'KR', 'XE', 'RN']
    # The upper limit for distances analysed for bonds will be assigned for
    # a given system (to save time). We take set('elements') and then find
    # the largest R(cov) in the system and set the max_dist as a double
    # of it plus the 150% tolerance (tol).
    set_of_elements = set(system['elements'])
    max_r_cov = max([
        pw.utilities.atomic_covalent_radius[i.upper()] for i in set_of_elements])
    max_dist = 2 * max_r_cov + tol
    # We continue untill all items in the list have been analysed and popped.
    count = 0
    already_out = []
    already_done = []
    II = ['C', 'C', 3.212, 4.085, 2.924]
    II2 = ['C', 'C', 1.703, 3.455, 8.276]
    II3 = ['C', 'C', 1.779, 3.22, 19.562]
    II4 = ['C', 'C', 1.629, 7.345, 1.516]
    II5 = ['C', 'C', 6.107, 9.313, 1.427]
    II6 = ['C', 'C', 5.216, 9.914, 5.734]
    II7 = ['C', 'C', -0.214, 9.925, 6.227]
    II8 = ['C', 'C', 0.205, 12.381, 0.758]
    II9 = ['C', 'C', 2.332, 14.01, 3.762]
    II10 = ['C', 'C', 9.328, 5.318, 2.113]
    II11 = ['C', 'C', 3.363, 6.383, 10.211]
    II12 = ['O', 'O', 7.69, 0.789, 8.018]
    II13 = ['C', 'C', 9.595, 4.169, 10.536]
    II14 = ['C', 'C', 9.507, 17.391, 2.006]
    II15 = ['C', 'C', 0.799, 16.503, 10.808]
    II16 = ['C', 'C', 9.839, 18.126, 10.293]
    II17 = ['O', 'O', 7.522, 23.228, 0.141]
    II18 = ['C', 'C', 7.006, 23.942, 5.255]
    II19 = ['C', 'C', 0.137, 25.585, 6.551]
    II20 = ['C', 'C', 0.449, 26.338, 0.515]
    II21 = ['C', 'C', 2.69, 29.263, 1.764]
    while atom_list:
        print('i', II in atom_list, II in satom_list)
        # print('i2', II2 in atom_list, II2 in satom_list)
        # print('i3', II3 in atom_list, II3 in satom_list)
        # print('i4', II4 in atom_list, II4 in satom_list)
        # print('i5', II5 in atom_list, II5 in satom_list)
        # print('i6', II6 in atom_list, II6 in satom_list)
        # print('i7', II7 in atom_list, II7 in satom_list)
        # print('i8', II8 in atom_list, II8 in satom_list)
        # print('i9', II9 in atom_list, II9 in satom_list)
        # print('i10', II10 in atom_list, II10 in satom_list)
        # print('i11', II11 in atom_list, II11 in satom_list)
        # print('i12', II12 in atom_list, II12 in satom_list)
        # print('i13', II13 in atom_list, II13 in satom_list)
        # print('i14', II14 in atom_list, II14 in satom_list)
        # print('i15', II15 in atom_list, II15 in satom_list)
        # print('i16', II16 in atom_list, II16 in satom_list)
        # print('i17', II17 in atom_list, II17 in satom_list)
        # print('i18', II18 in atom_list, II18 in satom_list)
        # print('i19', II19 in atom_list, II19 in satom_list)
        # print('i20', II20 in atom_list, II20 in satom_list)
        # print('i21', II21 in atom_list, II21 in satom_list)
        input()
        inside_atoms_heavy = [
            i for i in atom_list if i[0].upper() not in exceptions
        ]
        if inside_atoms_heavy:
            # Now we create an array of atom coordinates. It does seem
            # somehow counter-intuitive as this is what we started with
            # and made it into a list. But, in my opinion it's the only
            # way to do it. It's hard to control and delete items in two
            # separate arrays that we started with and we don't want
            # atoms already assigned in our array for distance matrix.
            inside_atoms_coord_heavy = pw.utilities.decompose_atom_list(inside_atoms_heavy)[
                1 + adj]
            dist_matrix = pw.utilities.euclidean_distances(inside_atoms_coord_heavy,
                                                           pseudo_origin.reshape(1, -1))
            atom_index_x, _ = np.unravel_index(dist_matrix.argmin(),
                                               dist_matrix.shape)
            # Added this so that lone atoms (even if heavy) close to the
            # periodic boundary are not analysed, as they surely have matching
            # symmetry equivalence that bind to a bigger atom cluster inside
            # the unit_cell.
            potential_starting_point = inside_atoms_heavy[atom_index_x]
            pot_arr = np.array(potential_starting_point[1 + adj:])
            dist_matrix = pw.utilities.euclidean_distances(
                atom_coor, pot_arr.reshape(1, -1)
            )
            idx = (dist_matrix > 0.1) * (dist_matrix < max_dist)
            if len(idx) < 1:
                print('2')
                pass
            else:
                working_list = [potential_starting_point]
                print('1')
            print(working_list)
            if potential_starting_point in already_done:
                print(potential_starting_point, 'done')
                input()
            already_done.append(potential_starting_point)
        else:
            # Safety check.
            break
        final_molecule = []
        while working_list:
            # print(working_list)
            # input()
            working_list_temp = []
            try:
                atom_coor = pw.utilities.decompose_atom_list(atom_list)[1 + adj]
            except _FunctionError:
                atom_coor = None
            for i in working_list:
                if i[0].upper() not in exceptions:
                    # It's of GREATEST importance that the i_arr variable
                    # is assigned here before entering the atom_coor loop.!
                    # Otherwise it will not be re-asigned when the satom_list
                    # still iterates, but the atom_list is already empty...
                    i_arr = np.array(i[1 + adj:])
                    if atom_coor is not None:
                        dist_matrix = pw.utilities.euclidean_distances(
                            atom_coor, i_arr.reshape(1, -1)
                        )
                        idx = (dist_matrix > 0.1) * (dist_matrix < max_dist)
                        neighbours_indexes = np.where(idx)[0]
                        for j in neighbours_indexes:
                            j_arr = np.array(atom_coor[j])
                            r_i_j = pw.utilities.distance(i_arr, j_arr)
                            r_cov_i_j = pw.utilities.atomic_covalent_radius[
                                i[0].upper()] + pw.utilities.atomic_covalent_radius[
                                    atom_list[j][0].upper()]
                            if r_cov_i_j - tol < r_i_j < r_cov_i_j + tol:
                                working_list_temp.append(atom_list[j])
                    if rebuild is not None:
                        sdist_matrix = pw.utilities.euclidean_distances(
                            satom_coor, i_arr.reshape(1, -1))
                        sidx = (sdist_matrix > 0.1) * (sdist_matrix < max_dist)
                        sneighbours_indexes = np.where(sidx)[0]
                        for j in sneighbours_indexes:
                            if satom_list[j] in atom_list:
                                pass
                            else:
                                j_arr = np.array(satom_coor[j])
                                r_i_j = pw.utilities.distance(i_arr, j_arr)
                                r_cov_i_j = pw.utilities.atomic_covalent_radius[
                                    i[0].upper()
                                ] + pw.utilities.atomic_covalent_radius[
                                    satom_list[j][0].upper()]
                                if r_cov_i_j - tol < r_i_j < r_cov_i_j + tol:
                                    working_list_temp.append(satom_list[j])
                    final_molecule.append(i)
                else:
                    final_molecule.append(i)
            print('ii', II in working_list)
            input()
            for i in working_list:
                # print(i in atom_list)
                if II in [i]:
                    print(i, II)
                    print(i in atom_list)
                    print('--')
                try:
                    atom_list.remove(i)
                    already_out.append(i)
                    if II in [i]:
                        print('removal')
                        print(i in atom_list)
                        while i in atom_list:
                            atom_list.remove(i)
                            print(i in atom_list)
                        sys.exit()
                except ValueError:
                    if II in [i]:
                        print('error!!')
                        print(i in atom_list)
                    # print(i in atom_list)
                    # print('s', i in already_out)
                    # print('---')
                    # print(i)
                    # # print(atom_list)
                    # # print(atom_coor)
                    # print('---')
                    if i not in atom_list:
                        pass
                    else:
                        sys.exit()
            # print('22', len(atom_list))
            # print(atom_list[:3])
            # input()
            # We empty the working list as all the items were analysed
            # and moved to the final_molecule list.
            working_list = []
            # We make sure there are no duplicates in the working_list_temp.
            working_list_temp = pw.utilities.unique(working_list_temp)
            # Now we move the entries from the temporary working list
            # to the working list for looping analysys.
            for i in working_list_temp:
                # We make sure that only new and unassigned atoms are
                # being transfered.
                if i not in final_molecule:
                    working_list.append(i)
        final_molecule_dict = {}
        final_molecule_dict['elements'] = np.array(
            [x[0] for x in final_molecule], dtype='str')
        final_molecule_dict['coordinates'] = np.array(
            [[*xyz[1 + adj:]] for xyz in final_molecule])
        if adj == 1:
            final_molecule_dict['atom_ids'] = np.array(
                [x[1] for x in final_molecule], dtype='str')

        # B = final_molecule_dict
        # D = ase.Atoms()
        # for i, j in enumerate(B['elements']):
        #     D.append(ase.Atom(position=B['coordinates'][i], symbol=j))
        # D.write(f'test_molecule_{count}.pdb')
        count += 1
        print(len(final_molecule_dict['atom_ids']))

        # In general we always want the molecule so the initial bool_ is True.
        bool_ = True
        # But, for periodic only if the molecule is in the initial unit cell.
        if rebuild is not None:
            com = pw.utilities.center_of_mass(final_molecule_dict['elements'],
                                              final_molecule_dict['coordinates'])
            com_frac = pw.utilities.fractional_from_cartesian(com, matrix)[0]
            # If we don't round the numerical errors will come up.
            com_frac_round = np.around(com_frac, decimals=8)
            bool_ = np.all(np.logical_and(com_frac_round >= boundary[0],
                                          com_frac_round < boundary[1]),
                           axis=0)
        bools.append(bool(bool_))
        if bool(bool_) is True:
            molecules.append(final_molecule_dict)
            # print(f'no molecules: {len(molecules)}')
    return molecules


def main():
    if (not len(sys.argv) == 1):
        print("""
    Usage: test_acdmnpbug.py
        """)
        sys.exit()

    file = 'ACDMNP_extracted.pdb'
    # file = 'TEST.pdb'
    # file = 'FOQTEM_extracted.pdb'

    logging.info(f'Run with modified functions')
    rbs = mod_mod(file=file)
    logging.info(f'{len(rbs)} molecules found!')
    sys.exit()

    logging.info(f'Run with pywindow functions')
    rbs = pywindow_f.modularize(file=file)
    Mol = rbs.molecules
    count = 0
    for molec in Mol:
        mol = Mol[molec]
        analysis = mol.full_analysis()
        if analysis is None:
            continue
        pdo = analysis['pore_diameter_opt']['diameter']
        if analysis['windows']['diameters'] is not None:
            nwind = len(analysis['windows']['diameters'])
        else:
            nwind = 0
        if pdo > 0.0 and nwind >= 2:
            # Mol[molec].dump_molecule(
            #     'testing' + "_MP_{0}_coms.pdb".format(molec),
            #     include_coms=True,
            #     override=True)
            # Mol[molec].dump_molecule(
            #     'testing' + "_MP_{0}.pdb".format(molec),
            #     include_coms=False,
            #     override=True)
            count += 1
    if count > 2:
        logging.info(f'Count = {count} -- too many molecules found.')

    for molec in Mol:
        mol = Mol[molec]
        logging.info(f'{mol.no_of_atoms}')


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='')
    main()
