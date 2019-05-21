#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to test pywindow results.

Author: Andrew Tarzia

Date Created: 20 May 2019

"""

import logging
import sys
import os
import numpy as np
import glob
import scipy.spatial.distance as scpy_dist
sys.path.insert(0, '/home/atarzia/thesource/')
import pywindow_f
import IO_tools


def main():
    if (not len(sys.argv) == 1):
        print("""
    Usage: pywindow_test.py
        """)
        sys.exit()
    else:
        pass

    CIFs = glob.glob('*cif')
    output_file = 'pywindow_test.out'
    output_str = ''
    # iterate over CIFs
    for cif in CIFs:
        file_prefix = cif.replace('.cif', '')
        output_str += cif + ':\n'
        pdb = IO_tools.convert_CIF_2_PDB(cif, wstruct=False)
        logging.info(f'> doing {pdb}')
        # modularize
        RB_s = pywindow_f.modularize(file=pdb)
        logging.info(f'> modularized')
        # run pywindow on each molecule
        for mol in RB_s.molecules:
            logging.info(f'> doing {mol}')
            if os.path.isfile(file_prefix + "_{0}.pdb".format(mol)):
                logging.info(f'> already done {mol}')
                continue
            Mol = RB_s.molecules[mol]
            analysis = Mol.full_analysis()
            print(analysis)
            # compare pore_diameter and pore_diameter_opt
            PD = analysis['pore_diameter']['diameter']
            PD_opt = analysis['pore_diameter_opt']['diameter']
            output_str += str(PD) + ','
            output_str += str(PD_opt) + ','
            logging.info(f'> PD {PD}, PD_opt {PD_opt}')
            # compare pore_volume and pore_volume_opt
            PV = analysis['pore_volume']
            PV_opt = analysis['pore_volume_opt']
            output_str += str(PV) + ','
            output_str += str(PV_opt) + ','
            logging.info(f'> PV {PV}, PV_opt {PV_opt}')
            # check min distance from cage COM to a cage atom
            COM = np.array([analysis['centre_of_mass']])
            coords = Mol.coordinates
            distances = scpy_dist.cdist(XA=COM, XB=coords, metric='euclidean')
            minD_C = min(distances[0])
            output_str += str(minD_C) + ','
            # check min distance from pore OPT COM to a cage atom
            PCOM = np.array([analysis['pore_diameter_opt']['centre_of_mass']])
            Pdistances = scpy_dist.cdist(XA=PCOM, XB=coords, metric='euclidean')
            minD_P = min(Pdistances[0])
            output_str += str(minD_P) + ','
            logging.info(f'> minD_C {minD_C}, minD_P {minD_P}')
            # output window number
            if analysis['windows']['diameters'] is not None:
                WN = len(analysis['windows']['diameters'])
            else:
                WN = 0
            output_str += str(WN) + ','
            logging.info(f'> WN {WN}')
            # wait for input - and visualize the structure
            Mol.dump_molecule(
                file_prefix + "_{0}.pdb".format(mol),
                include_coms=True,
                override=True)
            if PD_opt > 0:
                logging.info(f'> doing {pdb}')
                logging.info(f'> doing {mol}')
                logging.info(f'> PD {PD}, PD_opt {PD_opt}')
                logging.info(f'> PV {PV}, PV_opt {PV_opt}')
                logging.info(f'> minD_C {minD_C}, minD_P {minD_P}')
                logging.info(f'> WN {WN}')
                input()
            # INP = input('do these results make sense')
            INP = 'y'
            output_str += str(INP) + '\n'

    # output to file
    open(output_file, 'w').write(output_str)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='')
    main()
