#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to analyse crystal structures from PDB files.

Author: Andrew Tarzia

Date Created: 11 Nov 2020

"""


import sys

import stk

from utilities import (
    read_lib,
    get_plottables,
    calculate_ligand_SE,
    calculate_ligand_planarities,
    calculate_metal_ligand_distance,
)
from cage_analysis import write_xray_csv
from xtalcage import XtalCage


def main():
    first_line = (
        'Usage: analyse_crystal_structures.py '
        'complex_lib_file cage_set_lib_file ligand_directory '
        'cage_directory expt_lib_file'
    )
    if (not len(sys.argv) == 6):
        print(f"""
{first_line}

    complex_lib_file : (str)
        File containing complex information (XXXXX).

    cage_set_lib_file : (str)
        File containing cage information (XXXXX).

    ligand_directory : (str)
        Directory with required ligand structures.

    cage_directory : (str)
        Directory with required cage structures.

    expt_lib_file : (str)
        File containing experimental symmetry  information (XXXXX).

    """)
        sys.exit()
    else:
        complex_lib_file = sys.argv[1]
        cage_set_lib_file = sys.argv[2]
        ligand_directory = sys.argv[3]
        cage_directory = sys.argv[4]
        expt_lib_file = sys.argv[5]

    cage_set_lib = read_lib(cage_set_lib_file)
    complexes = read_lib(complex_lib_file)
    expt_data = read_lib(expt_lib_file)

    # List of the xtal structures and their corresponding names.
    xtals = {}
    for expt in expt_data:
        xtals[expt_data[expt]['xtal_struct_name']] = {
            'cage_set': expt,
            'symmetry_name': expt_data[expt]['symmetry'],
            'ligand_name': expt_data[expt]['ligand_name'],
            'complexes': tuple(expt_data[expt]['complexes']),
        }

    xtal_cage_data = {}
    comp_cage_data = {}
    for xtal in xtals:
        print(f'---- doing: {xtal}')
        pdb_file = f'{xtal}.pdb'
        cage_data = {}
        xtal_cage = XtalCage(
            name=xtal,
            pdb_file=pdb_file,
            complex_dicts=[
                complexes[i] for i in xtals[xtal]['complexes']
            ],
            cage_set_dict=cage_set_lib[xtals[xtal]['cage_set']]
        )
        org_ligs, smiles_keys = xtal_cage.get_organic_linkers()
        xtal_cage.write_metal_atom_structure()
        xtal_cage.collect_lowest_energy_conformer_file(
            cage_directory=cage_directory,
            n_atoms=[org_ligs[i].get_num_atoms() for i in org_ligs][0],
            cage_set=xtals[xtal]['cage_set'],
            ligand_name=xtals[xtal]['ligand_name'],
            ligand_directory=ligand_directory,
        )

        # Face-based analysis.
        m_structure = stk.BuildingBlock.init_from_file(
            f'{xtal_cage.name}_M.mol'
        )
        cage_data['m_cube_shape'] = (
             xtal_cage.get_m_shape(m_structure)
         )
        xtal_cage.define_faces(m_structure)
        cage_data['maxfacemetalpd'] = xtal_cage.get_max_face_metal_PD(
            m_structure
        )
        cage_data['maxintangledev'] = (
            xtal_cage.get_max_face_interior_angle_dev(m_structure)
        )
        cage_data['maxdifffaceaniso'] = (
            xtal_cage.get_max_face_anisotropy(m_structure)
        )

        # Full cage analysis.
        cage_data['porediam'] = xtal_cage.get_pore_size()
        cage_data['ML_lengths'] = calculate_metal_ligand_distance(
            mol=xtal_cage.stk_mol,
            metal_atomic_number=30,
            ligand_atomic_number=7,
        )
        cage_data['maxMLlength'] = max(cage_data['ML_lengths'])
        cage_data['octop'] = xtal_cage.get_min_order_value()

        # Ligand analysis.
        cage_data['core_planarities'] = calculate_ligand_planarities(
            org_ligs=org_ligs
        )
        cage_data['imine_torsions'] = (
            xtal_cage.calculate_abs_imine_torsions(org_ligs)
        )
        cage_data['strain_energies'] = calculate_ligand_SE(
            org_ligs=org_ligs,
            smiles_keys=smiles_keys,
            output_json=f'{xtal_cage.name}_lse.json',
            file_prefix=f'{xtal_cage.name}_sg'
        )
        cage_data['lsesum'] = sum([
            cage_data['strain_energies'][i]
            for i in cage_data['strain_energies']
        ])
        cage_data['minitors'] = min([
            j for i in cage_data['imine_torsions']
            for j in cage_data['imine_torsions'][i]
        ])
        cage_data['maxcrplan'] = max([
            cage_data['core_planarities'][i]
            for i in cage_data['core_planarities']
        ])
        xtal_cage_data[xtal] = cage_data
        comp_cage_data[xtal] = xtal_cage.get_cage_set_measures(
            cage_directory=cage_directory,
            cage_set=xtals[xtal]['cage_set'],
        )

        plottables = get_plottables(
            measures=comp_cage_data[xtal],
            name=xtal_cage.name
        )

        for p in plottables:
            p_dict = plottables[p]
            if p in ['formatione']:
                continue
            if p in ['lsesum']:
                xtal_data = cage_data[p] - min([
                    i for i in comp_cage_data[xtal][p].values()
                    if i is not None
                ])
            else:
                xtal_data = cage_data[p]

            expt_name = (
                f"C_{xtals[xtal]['cage_set']}_"
                f"{xtals[xtal]['symmetry_name']}"
            )
            xtal_cage.plot_Y(
                data=p_dict['data'],
                xtal_data=xtal_data,
                ylabel=p_dict['ylabel'],
                ylim=(None, None),
                filename=p_dict['filename'],
                expt_name=expt_name,
            )

    write_xray_csv(xtal_cage_data)


if __name__ == "__main__":
    main()
