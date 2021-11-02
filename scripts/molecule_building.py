#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Modules/functions for building molecules.

Author: Andrew Tarzia

Date Created: 23 Jan 2020

"""

import stk
import stko

from functional_groups import AromaticCNCFactory
import env_set


def build_oct_lambda(metal, ligand):

    complex = stk.ConstructedMolecule(
        topology_graph=stk.metal_complex.OctahedralLambda(
            metals={metal: 0},
            ligands={ligand: (0, 1, 2)},
        )
    )

    return complex


def build_oct_delta(metal, ligand):

    complex = stk.ConstructedMolecule(
        topology_graph=stk.metal_complex.OctahedralDelta(
            metals={metal: 0},
            ligands={ligand: (0, 1, 2)},
        )
    )

    return complex


def build_porphyrin(metal, ligand):

    complex = stk.ConstructedMolecule(
        topology_graph=stk.metal_complex.Porphyrin(
            metals={metal: 0},
            ligands={ligand: 0},
        )
    )

    return complex


def build_pw(metal, ligand):

    complex = stk.ConstructedMolecule(
        topology_graph=stk.metal_complex.Paddlewheel(
            metals={metal: (0, 1)},
            ligands={ligand: (0, 1, 2, 3)},
        )
    )

    return complex


def build_sqpl(metal, ligand):

    complex = stk.ConstructedMolecule(
        topology_graph=stk.metal_complex.SquarePlanar(
            metals={metal: 0},
            ligands={ligand: (0, 1, 2, 3)},
        )
    )

    return complex


def available_topologies(string):
    """
    Get function to build desired topology.

    """

    topologies = {
        'oct_lambda': build_oct_lambda,
        'oct_delta': build_oct_delta,
        'porphyrin': build_porphyrin,
        'pw': build_pw,
        'sqpl_mono': build_sqpl,
    }

    try:
        return topologies[string]
    except KeyError:
        raise KeyError(f'{string} not in {topologies.keys()}')


def custom_fg_factories(string):
    """
    Get factory of desired FG.

    """

    factories = {
        'pyridine_N_metal': AromaticCNCFactory(),
        'CO_metal': stk.SmartsFunctionalGroupFactory(
            smarts='[#6]~[#8X1]',
            bonders=(1, ),
            deleters=(),
        ),
        'COH_metal': stk.SmartsFunctionalGroupFactory(
            smarts='[#6]~[#8]~[#1]',
            bonders=(1, ),
            deleters=(2, ),
        ),
        'CNC_metal': stk.SmartsFunctionalGroupFactory(
            smarts='[#6]~[#7X2]~[#6]',
            bonders=(1, ),
            deleters=(),
        ),
        'CNBr_metal': stk.SmartsFunctionalGroupFactory(
            smarts='[#6]~[#7X2]~[#35]',
            bonders=(1, ),
            deleters=(),
        ),
    }

    try:
        return factories[string]
    except KeyError:
        raise KeyError(f'{string} not in {factories.keys()}')


def metal_FFs(CN):
    """
    Define metal FF names for UFF4MOF.

    Key = Atomic number
    Value = UFF4MOF type
    CN = coordination number of metal.

    """

    # Default settings.
    dicts = {
        26: 'Fe4+2',
        27: 'Co4+2',
        28: 'Ni4+2',  # No alternative available.
        30: 'Zn4+2',  # No alternative available for 90 degrees.
        42: 'Mo4f2',
        45: 'Rh6+3',  # No alternative available.
        46: 'Pd4+2',
        48: 'Cd4f2',  # No alternative available.
        78: 'Pt4+2',
    }

    if CN == 4:
        pass
    elif CN == 6:
        dicts[26] = 'Fe6+2'
        dicts[27] = 'Co6+2'

    return dicts


def optimize_SCA_complex(complex, name, dict, metal_FFs):
    """
    Optimize a sub-component self assmebly complex.

    """

    print(f'doing UFF4MOF optimisation for {name}')
    gulp_opt = stko.GulpUFFOptimizer(
        gulp_path=env_set.gulp_path(),
        metal_FF=metal_FFs,
        output_dir=f'{name}_uff1'
    )
    gulp_opt.assign_FF(complex)
    complex = gulp_opt.optimize(mol=complex)
    complex.write(f'{name}_uff1.mol')

    print(f'doing xTB optimisation for {name}')
    xtb_opt = stko.XTB(
        xtb_path=env_set.xtb_path(),
        output_dir=f'{name}_xtb',
        gfn_version=2,
        num_cores=6,
        opt_level='tight',
        charge=dict['total_charge'],
        num_unpaired_electrons=dict['unpaired_e'],
        max_runs=1,
        calculate_hessian=False,
        unlimited_memory=True
    )
    complex = xtb_opt.optimize(mol=complex)
    complex.write(f'{name}_opt.mol')

    return complex
