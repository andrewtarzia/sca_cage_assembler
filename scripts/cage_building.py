#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Modules defining and building the Cage class and subclasses.

Author: Andrew Tarzia

Date Created: 03 Mar 2020

"""

from itertools import product
from os.path import exists, join
import json
import pywindow as pw
import os

import stk

import atools
from molecule_building import metal_FFs


def available_topologies(string):
    """
    Get stk function of desired topology.

    """

    topologies = {
        'm4l4spacer': stk.cage.M4L4_Oct_Spacer(),
        'm8l6face': stk.cage.M8L6_Oct_Face(),
        'm6l2l3': stk.cage.M6L2L3_Oct()
    }

    try:
        return string, topologies[string]
    except KeyError:
        raise KeyError(f'{string} not in {topologies.keys()}')


class Cage:
    """
    Generic class that builds and analyses stk.ConstructuedMolecules.

    """

    def __init__(
        self,
        name,
        bbs,
        topology,
        bb_vertices,
        charge,
        free_electron_options
    ):

        self.name = name
        self.bbs = bbs
        self.topology = topology
        self.bb_vertices = bb_vertices
        self.unopt_file = f'{self.name}_unopt'
        self.crush_file = f'{self.name}_cru'
        self.uff4mof_file = f'{self.name}_uff'
        self.uffMD_file = f'{self.name}_prextb'
        self.opt_file = f'{self.name}_optc'
        self.pw_file = f'{self.name}_pw'
        self.op_file = f'{self.name}_OP'
        self.charge = charge
        self.free_electron_options = free_electron_options
        print(self.charge, self.free_electron_options)

    def build(self):
        print(f'....building {self.name}')
        cage = stk.ConstructedMolecule(
            building_blocks=self.bbs,
            topology_graph=self.topology,
            building_block_vertices=self.bb_vertices
        )
        cage.write(f'{self.unopt_file}.mol')
        cage.dump(f'{self.unopt_file}.json')
        self.cage = cage

    def optimize(self, free_e):
        custom_metal_FFs = metal_FFs(CN=6)

        # Skip if _opt.json exists.
        if exists(f'{self.opt_file}.json'):
            self.cage.update_from_file(f'{self.opt_file}.mol')
            return
        print(f'....optimizing {self.name}')

        # Run if crush output does not exist.
        if not exists(f'{self.crush_file}.mol'):
            self.cage = atools.MOC_collapse(
                self.cage,
                self.name,
                step_size=0.05,
                distance_cut=2.0
            )
            self.cage.write(f'{self.crush_file}.mol')
        else:
            self.cage.update_from_file(f'{self.crush_file}.mol')

        # Run if uff4mof opt output does not exist.
        if not exists(f'{self.uff4mof_file}.mol'):
            self.cage = atools.MOC_uff_opt(
                self.cage,
                self.name,
                metal_FFs=custom_metal_FFs
            )
            self.cage.write(f'{self.uff4mof_file}.mol')
        else:
            self.cage.update_from_file(f'{self.uff4mof_file}.mol')

        # Run if uff4mof MD output does not exist.
        if not exists(f'{self.uffMD_file}.mol'):
            self.cage = atools.MOC_MD_opt(
                self.cage,
                self.name,
                integrator='leapfrog verlet',
                temperature='700',
                N=10,
                timestep='0.5',
                equib='0.1',
                production='2',
                metal_FFs=custom_metal_FFs,
                opt_conf=False,
                save_conf=False
            )
            self.cage.write(f'{self.uffMD_file}.mol')
        else:
            self.cage.update_from_file(f'{self.uffMD_file}.mol')

        atools.MOC_xtb_opt(
            self.cage,
            self.name,
            nc=6,
            free_e=free_e,
            charge=self.charge,
            opt_level='normal',
            etemp=300,
            # solvent=('dmso', 'verytight')
        )
        self.cage.write(f'{self.opt_file}.mol')
        self.cage.dump(f'{self.opt_file}.json')

    def get_energy(self):
        print(f'....getting energy of {self.name}')
        raise NotImplementedError()

    def compare_UHF_values(self):
        print(f'....comparing UHF of {self.name}')
        raise NotImplementedError()

    def analyze_cage_geometry(self):
        """
        Analyse cage geometry using order parameters.

        """

        # Check if output file exists.
        if not exists(f'{self.op_file}.json'):
            print(f'....analyzing geometry of {self.name}')

            # Get atomic numbers of all present metals.
            pres_atm_no = list(set([
                i.atomic_number for i in self.cage.atoms
                if i.atomic_number in metal_FFs(CN=4).keys()
            ]))

            # Get OPs for each metal independantly.
            # Save to dict with atom id and atom type.
            op_res = {}
            for metal in pres_atm_no:
                op_set = atools.get_order_values(
                    mol=self.cage,
                    metal=metal,
                    per_site=True
                )
                op_res[metal] = op_set

            # Save to output files.
            with open(f'{self.op_file}.json', 'w') as f:
                json.dump(op_res, f)

        # Get data.
        with open(f'{self.op_file}.json', 'r') as f:
            self.op_data = json.load(f)

    def analyze_cage_porosity(self):
        """
        Analyse cage porosity with pywindow.

        """

        # Check if output file exists.
        if not exists(f'{self.pw_file}.json'):
            print(f'....analyzing porosity of {self.name}')
            # Load cage into pywindow.
            self.cage.write('temp.xyz')
            pw_cage = pw.MolecularSystem.load_file(
                'temp.xyz'
            )
            pw_cage_mol = pw_cage.system_to_molecule()
            os.system('rm temp.xyz')

            # Calculate pore size.
            pw_cage_mol.calculate_pore_diameter_opt()
            pw_cage_mol.calculate_pore_volume_opt()

            # Save files.
            pw_cage_mol.dump_properties_json(f'{self.pw_file}.json')
            pw_cage_mol.dump_molecule(
                f'{self.pw_file}.pdb',
                include_coms=True
            )

        # Get data.
        with open(f'{self.pw_file}.json', 'r') as f:
            self.pw_data = json.load(f)

    def __str__(self):
        return (
            f'{self.__class__.__name__}'
            f'(name={self.name}, topology={self.topology})'
        )

    def __repr__(self):
        return str(self)


class HetPrism:
    """
    Class that builds and analyses stk.ConstructuedMolecules.

    Represents heteroleptic prism cages and all necessary homoleptic
    cages.

    """

    def __init__(
        self,
        name,
        prism_dict,
        complex_dicts,
        ligand_dir,
        complex_dir
    ):
        self.name = name
        self.properties_file = f'{self.name}_HP.json'
        self.prism_dict = prism_dict
        self.complex_dicts = complex_dicts
        self.cages_to_build = self.define_cages_to_build(
            ligand_dir,
            complex_dir
        )
        self.built_cage_properties = {}

    def load_properties(self):
        """
        Load class from JSON file.

        """
        with open(self.properties_file, 'r') as f:
            self.built_cage_properties = json.load(f)

    def dump_properties(self):
        """
        Dump class to JSON file.

        """
        with open(self.properties_file, 'w') as f:
            json.dump(self.built_cage_properties, f)

    def __str__(self):
        return (
            f'{self.__class__.__name__}'
            f'(name={self.name})\n'
            f'{self.prism_dict}'
        )

    def __repr__(self):
        return str(self)

    def get_ratios(self, n_metals):
        rng = range(0, n_metals+1)
        rats = []
        for i in product(rng, rng):
            if i[0]+i[1] == n_metals:
                rats.append(i)
        return rats

    def load_complex(self, complex_name, complex_dir):
        unopt_c = stk.ConstructedMolecule.load(
            join(complex_dir, f'{complex_name}_opt.json')
        )
        complex = stk.BuildingBlock.init_from_molecule(
            unopt_c,
            functional_groups=['bromine']
        )

        return complex

    def load_ligand(self, ligand_name, ligand_dir):
        unopt_l = stk.Molecule.load(
            join(ligand_dir, f'{ligand_name}_opt.json')
        )
        ligand = stk.BuildingBlock.init_from_molecule(
            unopt_l,
            functional_groups=['bromine']
        )

        with open(join(ligand_dir, 'ligand_data.json'), 'r') as f:
            lig_prop = json.load(f)
        properties = lig_prop[ligand_name]

        return ligand, properties

    def define_cages_to_build(self, ligand_dir, complex_dir):
        """
        Defines the name and objects of all cages to build.

        """

        cages_to_build = []

        # Get Delta and Lambda complex.
        print(self.complex_dicts)
        D_complex_name = [
            i for i in self.complex_dicts if 'del' in i
        ][0]
        L_complex_name = [
            i for i in self.complex_dicts if 'lam' in i
        ][0]
        print(D_complex_name, L_complex_name)
        L_complex = self.load_complex(
            complex_name=L_complex_name,
            complex_dir=complex_dir
        )
        D_complex = self.load_complex(
            complex_name=D_complex_name,
            complex_dir=complex_dir
        )

        D_charge = self.complex_dicts[D_complex_name]['total_charge']
        L_charge = self.complex_dicts[L_complex_name]['total_charge']
        D_free_e = self.complex_dicts[D_complex_name][
            'unpaired_e'
        ].strip(')()').split(',')
        L_free_e = self.complex_dicts[L_complex_name][
            'unpaired_e'
        ].strip(')()').split(',')
        print(D_charge, D_free_e, L_charge, L_free_e)

        # Get all linkers.
        tet_linker, tet_prop = self.load_ligand(
            ligand_name=self.prism_dict['tetratopic'],
            ligand_dir=ligand_dir
        )
        tri_linker, tri_prop = self.load_ligand(
            ligand_name=self.prism_dict['tritopic'],
            ligand_dir=ligand_dir
        )

        print(D_complex, L_complex, tet_linker, tri_linker)

        # Tetratopic homoleptic cages of all symmetries.
        tet_topo_name, tet_topo = available_topologies(
            string='m8l6face'
        )

        tet_n_metals = 8
        tet_ratios = self.get_ratios(tet_n_metals)
        for rat in tet_ratios:
            print(rat)
            new_name = (
                f"C_{self.prism_dict['corner_name']}_"
                f"{self.prism_dict['tetratopic']}_"
                f"{tet_topo_name}_"
                f"d{rat[0]}l{rat[1]}_"
                f"SYMM"
            )
            new_bbs = [D_complex, L_complex, tet_linker]
            new_bb_vertices = {
                D_complex: tet_topo.vertices[:rat[0]],
                L_complex: tet_topo.vertices[rat[0]:rat[0]+rat[1]],
                tet_linker: tet_topo.vertices[tet_n_metals:]
            }
            # Merge linker and complex charges.
            complex_charge = rat[0]*int(D_charge)+rat[1]*int(L_charge)
            new_charge = int(tet_prop[0])*6 + complex_charge

            print(tet_prop[1])
            lig_free_e = [int(i)*6 for i in tet_prop[1]]
            compl_free_e = [
                int(i)*rat[0] + int(j)*rat[1]
                for i, j in zip(D_free_e, L_free_e)
            ]
            print(lig_free_e, compl_free_e)

            new_free_electron_options = []
            for opt in product(lig_free_e, compl_free_e):
                print(opt)
                new_free_electron_options.append(opt[0]+opt[1])

            print(new_charge, new_free_electron_options)
            new_cage = Cage(
                name=new_name,
                bbs=new_bbs,
                topology=tet_topo,
                bb_vertices=new_bb_vertices,
                charge=new_charge,
                free_electron_options=new_free_electron_options
            )
            cages_to_build.append(new_cage)
            print('NOT BUILDING ALL RATIOS CURRENTLY!!!!!!')
            break

        # Tritopic homoleptic cages of all symmetries.
        tri_topo_name, tri_topo = available_topologies(
            string='m4l4spacer'
        )

        tri_n_metals = 4
        tri_ratios = self.get_ratios(tri_n_metals)
        for rat in tri_ratios:
            print(rat)
            new_name = (
                f"C_{self.prism_dict['corner_name']}_"
                f"{self.prism_dict['tritopic']}_"
                f"{tri_topo_name}_"
                f"d{rat[0]}l{rat[1]}_"
                f"SYMM"
            )
            new_bbs = [D_complex, L_complex, tri_linker]
            new_bb_vertices = {
                D_complex: tet_topo.vertices[:rat[0]],
                L_complex: tet_topo.vertices[rat[0]:rat[0]+rat[1]],
                tri_linker: tet_topo.vertices[tri_n_metals:]
            }
            # Merge linker and complex charges.
            complex_charge = rat[0]*int(D_charge)+rat[1]*int(L_charge)
            new_charge = int(tri_prop[0])*4 + complex_charge

            print(tri_prop[1])
            lig_free_e = [int(i)*4 for i in tri_prop[1]]
            compl_free_e = [
                int(i)*rat[0] + int(j)*rat[1]
                for i, j in zip(D_free_e, L_free_e)
            ]
            print(lig_free_e, compl_free_e)

            new_free_electron_options = []
            for opt in product(lig_free_e, compl_free_e):
                print(opt)
                new_free_electron_options.append(opt[0]+opt[1])

            print(new_charge, new_free_electron_options)
            new_cage = Cage(
                name=new_name,
                bbs=new_bbs,
                topology=tri_topo,
                bb_vertices=new_bb_vertices,
                charge=new_charge,
                free_electron_options=new_free_electron_options
            )
            cages_to_build.append(new_cage)
            print('NOT BUILDING ALL RATIOS CURRENTLY!!!!!!')
            break

        # Prisms of all symmetries.
        pri_topo_name, pri_topo = available_topologies(
            string='m6l2l3'
        )

        pri_n_metals = 6
        pri_ratios = self.get_ratios(pri_n_metals)
        for rat in pri_ratios:
            new_name = (
                f"HeP_{self.prism_dict['corner_name']}_"
                f"{self.prism_dict['tritopic']}_"
                f"{self.prism_dict['tetratopic']}_"
                f"{pri_topo_name}_"
                f"d{rat[0]}l{rat[1]}_"
                f"SYMM"
            )
            new_bbs = [D_complex, L_complex, tri_linker, tet_linker]
            new_bb_vertices = {
                D_complex: tet_topo.vertices[:rat[0]],
                L_complex: tet_topo.vertices[rat[0]:rat[0]+rat[1]],
                tri_linker: tet_topo.vertices[
                    pri_n_metals:pri_n_metals+2
                ],
                tet_linker: tet_topo.vertices[pri_n_metals+2:]
            }
            # Merge linker and complex charges.
            complex_charge = rat[0]*int(D_charge)+rat[1]*int(L_charge)
            new_charge = int(tri_prop[0])*2
            new_charge += int(tet_prop[0])*3 + complex_charge

            print(tri_prop[1], tet_prop[1])
            tri_free_e = [int(i)*2 for i in tri_prop[1]]
            tet_free_e = [int(i)*3 for i in tet_prop[1]]
            compl_free_e = [
                int(i)*rat[0] + int(j)*rat[1]
                for i, j in zip(D_free_e, L_free_e)
            ]
            print(tri_free_e, tet_free_e, compl_free_e)

            new_free_electron_options = []
            for opt in product(tri_free_e, tet_free_e, compl_free_e):
                print(opt)
                new_free_electron_options.append(opt[0]+opt[1]+opt[2])

            print(new_charge, new_free_electron_options)
            new_cage = Cage(
                name=new_name,
                bbs=new_bbs,
                topology=pri_topo,
                bb_vertices=new_bb_vertices,
                charge=new_charge,
                free_electron_options=new_free_electron_options
            )
            cages_to_build.append(new_cage)
            print('NOT BUILDING ALL RATIOS CURRENTLY!!!!!!')
            break

        return cages_to_build
