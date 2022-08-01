#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot Zn-Zn distributions in calc'd and xrd structures.

Author: Andrew Tarzia

"""

import stk
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
from utilities import read_lib
import scipy.spatial.distance


def get_pairwise_dists(stk_mol):
    dists = scipy.spatial.distance.pdist(
        stk_mol.get_position_matrix(),
    )
    return dists


def main():
    first_line = (
        'Usage: plot_znzn_distributions.py xray_structure_path expt_lib_file'
    )
    if (not len(sys.argv) == 3):
        print(f"""
{first_line}

    xray_structure_path : (str)
        ../xray_structures/analysis/

    expt_lib_file : (str)
        File containing experimental symmetry  information (XXXXX).

    """)
        sys.exit()
    else:
        xray_structure_path = sys.argv[1]
        expt_lib_file = sys.argv[2]

    expt_data = read_lib(expt_lib_file)

    _figure_path = 'figures'

    for xray_name in expt_data:
        cs = xray_name.split('-')[0]
        cs_data = expt_data[xray_name]
        symm = cs_data['symmetry']
        xname = cs_data['xtal_struct_name']
        xtal_file = os.path.join(
            xray_structure_path,
            f'{xname}_M.mol'
        )
        meta_file = f'{cs}_{symm}_optc_M.mol'
        calc_file = f'C_{cs}_{symm}_optc.mol'
        xtal_structure = stk.BuildingBlock.init_from_file(xtal_file)
        calc_structure = stk.BuildingBlock.init_from_file(calc_file)
        metal_atom_ids = [
            i.get_id() for i in calc_structure.get_atoms()
            if i.get_atomic_number() == 30
        ]

        # Write to mol file.
        calc_structure.write(meta_file, atom_ids=metal_atom_ids)
        meta_structure = stk.BuildingBlock.init_from_file(meta_file)
        print(meta_structure)
        print(xtal_structure)
        meta_pairs = get_pairwise_dists(meta_structure)
        print(meta_pairs)
        xtal_pairs = get_pairwise_dists(xtal_structure)
        print(xtal_pairs)

        xmin = 9
        xmax = 16.5
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.hist(
            meta_pairs,
            bins=np.arange(xmin, xmax, 0.2),
            density=False,
            color='k',
            alpha=1.0,
            histtype='step',
            lw=2,
            label='calculated'
        )
        ax.hist(
            xtal_pairs,
            bins=np.arange(xmin, xmax, 0.2),
            density=False,
            color='#e71989',
            alpha=1.0,
            histtype='step',
            lw=2,
            label='xray'
        )
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel(r'distance [$\mathrm{\AA}$]', fontsize=16)
        ax.set_ylabel('count', fontsize=16)
        fig.legend(fontsize=16)

        fig.tight_layout()
        fig.savefig(
            os.path.join(_figure_path, f"znzndist_{xray_name}.pdf"),
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()


if __name__ == "__main__":
    main()
