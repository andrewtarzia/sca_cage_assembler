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
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import json

import stk

import atools
import symmetries
from utilities import calculate_binding_AR
from cage_building import available_topologies
from cage import Cage


class CageSet:
    """
    Class that builds and analyses stk.ConstructuedMolecules.

    """

    def __init__(
        self,
        name,
        cage_set_dict,
        complex_dicts,
        ligand_dicts,
        ligand_dir,
        complex_dir
    ):
        self.name = name
        self.properties_file = f'{self.name}_CS.json'
        self.cage_set_dict = cage_set_dict
        self.complex_dicts = complex_dicts
        self.ligand_dicts = ligand_dicts
        self.cages_to_build = self.define_cages_to_build(
            ligand_dir,
            complex_dir
        )
        self.built_cage_properties = {}

    def define_cages_to_build(self):
        """
        Defines the name and objects of all cages to build.

        """

        raise NotImplementedError(
            f'Not implemented for {self.__class__}'
        )

    def load_properties(self):
        """
        Load class from JSON file.

        """

        if exists(self.properties_file):
            with open(self.properties_file, 'r') as f:
                self.built_cage_properties = json.load(f)

    def dump_properties(self):
        """
        Dump class to JSON file.

        """
        with open(self.properties_file, 'w') as f:
            json.dump(self.built_cage_properties, f, indent=4)

    def __str__(self):
        return (
            f'{self.__class__.__name__}'
            f'(name={self.name})\n'
            f'{self.cage_set_dict}'
        )

    def __repr__(self):
        return str(self)

    def _get_no_vertices(self, string):
        """
        Get the number of vertices for a given topology.

        """

        topologies = {
            'm4l4spacer': 8,
            'm8l6face': 14,
            'm6l2l3': 11,
        }

        try:
            return topologies[string]
        except KeyError:
            raise KeyError(f'{string} not in {topologies.keys()}')

    def _get_complex_info(self, complex_dir):

        D_complex_name = [
            i for i in self.complex_dicts if 'del' in i
        ][0]
        L_complex_name = [
            i for i in self.complex_dicts if 'lam' in i
        ][0]
        L_complex = self._load_complex(
            complex_name=L_complex_name,
            complex_dir=complex_dir
        )
        D_complex = self._load_complex(
            complex_name=D_complex_name,
            complex_dir=complex_dir
        )

        return D_complex_name, D_complex, L_complex_name, L_complex

    def _get_complex_properties(self, complex_name):
        charge = self.complex_dicts[complex_name]['total_charge']
        free_e = self.complex_dicts[complex_name]['unpaired_e']

        return charge, free_e

    def _get_ligand(self, type_name, ligand_dir):
        prop = self.ligand_dicts[self.cage_set_dict[type_name]]
        linker = self._load_ligand(
            ligand_name=self.cage_set_dict[type_name],
            ligand_dir=ligand_dir
        )

        return prop, linker

    def _get_rot_vertices(self, string):
        """
        Get the list of rotatable vertices for a given topology.

        Only ligand vertices are rotatable in this case.

        # TODO: Currently only defined for cube (90 deg). Add tri-face.

        """

        if string in ['m4l4spacer', 'm6l2l3']:
            raise NotImplementedError(
                'Currently only defined for cube (90 deg). Add tri.'
            )

        topologies = {
            'm4l4spacer': [4, 5, 6, 7],
            'm8l6face': [8, 9, 10, 11, 12, 13],
            'm6l2l3': [8, 9, 10]
        }

        try:
            return topologies[string]
        except KeyError:
            raise KeyError(f'{string} not in {topologies.keys()}')

    def _get_ratios(self, n_metals):
        rng = range(0, n_metals+1)
        rats = []
        for i in product(rng, rng):
            if i[0]+i[1] == n_metals:
                rats.append(i)
        return rats

    def _load_complex(self, complex_name, complex_dir):
        complex = stk.BuildingBlock.init_from_file(
            join(complex_dir, f'{complex_name}_opt.mol'),
            functional_groups=[stk.BromoFactory()]
        )

        return complex

    def _load_ligand(self, ligand_name, ligand_dir):
        ligand = stk.BuildingBlock.init_from_file(
            join(ligand_dir, f'{ligand_name}_opt.mol'),
            functional_groups=[stk.BromoFactory()]
        )

        return ligand

    def get_cage_symmetries(
        self,
        string,
        D_complex,
        L_complex,
        linkers,
    ):
        """
        Returns cage symmetries for a given topology.

        """

        if string == 'm4l4spacer':
            symm_list = {}
            linker = linkers[3]

            # Predefined list of symmetries.
            symm_c = symmetries.M4L4_Symmetry(
                D_complex=D_complex,
                L_complex=L_complex,
                linker=linker,
            )
            symm_list['base'] = symm_c.base()

        elif string == 'm8l6face':
            symm_list = {}
            linker = linkers[4]

            # Predefined list of symmetries.
            symm_c = symmetries.M8L6_Symmetry(
                D_complex=D_complex,
                L_complex=L_complex,
                linker=linker,
            )
            symm_list['o1'] = symm_c.o1()
            symm_list['th1'] = symm_c.th1()
            symm_list['th2'] = symm_c.th2()
            symm_list['t1'] = symm_c.t1()
            symm_list['s61'] = symm_c.s61()
            symm_list['s62'] = symm_c.s62()
            symm_list['d31'] = symm_c.d31()
            symm_list['d32'] = symm_c.d32()
            symm_list['c2v'] = symm_c.c2v()
            symm_list['c2h'] = symm_c.c2h()

        elif string == 'm6l2l3':
            symm_list = {}
            linker3 = linkers[3]
            linker4 = linkers[4]

            # Predefined list of symmetries.
            symm_c = symmetries.M6L2L3_Symmetry(
                D_complex=D_complex,
                L_complex=L_complex,
                linker3=linker3,
                linker4=linker4,
            )
            symm_list['base'] = symm_c.base()
        else:
            raise KeyError(f'{string} not in defined')

        print(f'{len(symm_list)} symmetries to build')
        return symm_list


class HoCube(CageSet):
    """
    Class that builds and analyses stk.ConstructuedMolecules.

    Represents homoleptic cube cages with all necessary symmetries
    and orientations.

    """

    def __init__(
        self,
        name,
        cage_set_dict,
        complex_dicts,
        ligand_dicts,
        ligand_dir,
        complex_dir
    ):

        super().__init__(
            name,
            cage_set_dict,
            complex_dicts,
            ligand_dicts,
            ligand_dir,
            complex_dir
        )

        # Get ligand aspect ratio.
        self.ligand_aspect_ratio = self._get_ligand_AR(ligand_dir)

    def _get_ligand_AR(self, ligand_dir):
        """
        Calculate ligand aspect ratio based on binder positions.

        """

        tet_linker = self._load_ligand(
            ligand_name=self.cage_set_dict['tetratopic'],
            ligand_dir=ligand_dir
        )
        ligand_AR = calculate_binding_AR(tet_linker)

        return ligand_AR

    def define_cages_to_build(self, ligand_dir, complex_dir):
        """
        Defines the name and objects of all cages to build.

        """

        # Get Delta and Lambda complexes.
        D_complex_name, D_complex, L_complex_name, L_complex = (
            self._get_complex_info(complex_dir=complex_dir)
        )

        D_charge, D_free_e = self._get_complex_properties(
            complex_name=D_complex_name
        )
        L_charge, L_free_e = self._get_complex_properties(
            complex_name=L_complex_name
        )

        # Get linker and dictionary.
        tet_prop, tet_linker = self._get_ligand(
            type_name='tetratopic',
            ligand_dir=ligand_dir
        )

        # Get topology function as object to be used in following list.
        # Homoleptic cage with tetratopic ligand.
        tet_topo_name = 'm8l6face'
        tet_topo_fn = available_topologies(string=tet_topo_name)

        symmetries_to_build = self.get_cage_symmetries(
            string=tet_topo_name,
            D_complex=D_complex,
            L_complex=L_complex,
            linkers={4: tet_linker},
        )

        for name_string in symmetries_to_build:
            base_name = (
                f"C_{self.cage_set_dict['corner_name']}_"
                f"{self.cage_set_dict['tetratopic']}"
            )
            new_name = f"{base_name}_{name_string}"
            building_blocks = (
                symmetries_to_build[name_string]['building_blocks']
            )
            vertex_alignments = (
                symmetries_to_build[name_string]['vertex_alignments']
            )
            rat = symmetries_to_build[name_string]['ratio']

            # Merge linker and complex charges.
            complex_charge = rat[0]*int(D_charge)
            complex_charge += rat[1]*int(L_charge)
            new_charge = tet_prop['net_charge']*6 + complex_charge

            lig_free_e = [
                int(i)*6 for i in tet_prop['total_unpaired_e']
            ]
            compl_free_e = [
                int(i)*rat[0] + int(j)*rat[1]
                for i, j in zip(D_free_e, L_free_e)
            ]
            new_free_electron_options = []
            for opt in product(lig_free_e, compl_free_e):
                new_free_electron_options.append(opt[0]+opt[1])

            new_cage = Cage(
                name=new_name,
                base_name=base_name,
                topology_fn=tet_topo_fn,
                building_blocks=building_blocks,
                vertex_alignments=vertex_alignments,
                topology_string=tet_topo_name,
                charge=new_charge,
                free_electron_options=new_free_electron_options,
                cage_set_dict=self.cage_set_dict
            )
            cages_to_build.append(new_cage)

        return cages_to_build

    def plot_Y(self, data, ylabel, filename, ylim=None):
        C = '#AFE074'
        M = 'o'

        fig, ax = plt.subplots(figsize=(8, 5))
        x_pos_list = []
        names_list = []
        for i, name in enumerate(data):
            X = i+2
            names_list.append(name.split('_')[-1])
            x_pos_list.append(X)
            ax.scatter(
                X,
                data[name],
                c=C,
                edgecolors='k',
                marker=M,
                alpha=1.0,
                s=180
            )
        # Set number of ticks for x-axis
        ax.tick_params(axis='both', which='major', labelsize=16)
        # ax.set_xlabel(r'pore volume [$\mathrm{\AA}^3$]', fontsize=16)
        ax.set_ylabel(ylabel, fontsize=16)
        ax.set_xlim(1, i+3)
        ax.set_ylim(ylim)
        ax.set_xticklabels(names_list)
        ax.set_xticks(x_pos_list)

        fig.tight_layout()
        fig.savefig(
            filename,
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()

    def plot_Y_C(
        self,
        data,
        ylabel,
        data_C,
        clabel,
        clim,
        filename,
        ylim=None,
    ):

        M = 'o'

        fig, ax = plt.subplots(figsize=(8, 5))
        # Define cmap.
        CMAP = {i: data_C[i]/clim[1] for i in data_C}
        cmap = {
            'mid_point': 0.5,
            'cmap': cm.Purples_r,
            'ticks': [0, .50, 1.00],
            'labels': [
                str(clim[0]),
                str((clim[1]-clim[0])/2),
                str(clim[1])
            ],
            'cmap_label': clabel,
        }

        cmp = atools.define_plot_cmap(
            fig, ax,
            mid_point=cmap['mid_point'],
            cmap=cmap['cmap'],
            ticks=cmap['ticks'],
            labels=cmap['labels'],
            cmap_label=cmap['cmap_label']
        )

        x_pos_list = []
        names_list = []
        xs = []
        ys = []
        cs = []
        for i, name in enumerate(data):
            X = i+2
            names_list.append(name.split('_')[-1])
            x_pos_list.append(X)
            xs.append(X)
            ys.append(data[name])
            cs.append(cmp(CMAP[name]))
        ax.scatter(
            xs,
            ys,
            c=cs,
            edgecolors='k',
            marker=M,
            alpha=1.0,
            s=180
        )

        # Set number of ticks for x-axis
        ax.tick_params(axis='both', which='major', labelsize=16)
        # ax.set_xlabel(r'pore volume [$\mathrm{\AA}^3$]', fontsize=16)
        ax.set_ylabel(ylabel, fontsize=16)
        ax.set_xlim(1, i+3)
        ax.set_ylim(ylim)
        ax.set_xticklabels(names_list)
        ax.set_xticks(x_pos_list)

        fig.tight_layout()
        fig.savefig(
            filename,
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()


class HetPrism(CageSet):
    """
    Class that builds and analyses stk.ConstructuedMolecules.

    Represents heteroleptic prism cages and all necessary homoleptic
    cages.

    """

    def __init__(
        self,
        name,
        cage_set_dict,
        complex_dicts,
        ligand_dicts,
        ligand_dir,
        complex_dir
    ):

        super().__init__(
            name,
            cage_set_dict,
            complex_dicts,
            ligand_dicts,
            ligand_dir,
            complex_dir
        )

    def define_cages_to_build(self, ligand_dir, complex_dir):
        """
        Defines the name and objects of all cages to build.

        """

        # Get Delta and Lambda complex.
        D_complex_name, D_complex, L_complex_name, L_complex = (
            self._get_complex_info(complex_dir=complex_dir)
        )

        D_charge, D_free_e = self._get_complex_properties(
            complex_name=D_complex_name
        )
        L_charge, L_free_e = self._get_complex_properties(
            complex_name=L_complex_name
        )

        # Get linker and dictionary.
        tet_prop, tet_linker = self._get_ligand(
            type_name='tetratopic',
            ligand_dir=ligand_dir
        )
        tri_prop, tri_linker = self._get_ligand(
            type_name='tritopic',
            ligand_dir=ligand_dir
        )

        # Get topology function as object to be used in following list.
        # Homoleptic cage with tetratopic ligand.
        tet_topo_name = 'm8l6face'
        tet_topo_fn = available_topologies(string=tet_topo_name)
        # Homoleptic cage with tritopic ligand.
        tri_topo_name = 'm4l4spacer'
        tri_topo_fn = available_topologies(string=tri_topo_name)
        # Heteroleptic cage with tetratopic + tritopic ligand.
        pri_topo_name = 'm6l2l3'
        pri_topo_fn = available_topologies(string=pri_topo_name)

        # Define symmetries of each topology.
        tet_symmetries_to_build = self.get_cage_symmetries(
            string=tet_topo_name,
            D_complex=D_complex,
            L_complex=L_complex,
            linkers={4: tet_linker},
        )
        tri_symmetries_to_build = self.get_cage_symmetries(
            string=tri_topo_name,
            D_complex=D_complex,
            L_complex=L_complex,
            linkers={3: tri_linker},
        )
        pri_symmetries_to_build = self.get_cage_symmetries(
            string=pri_topo_name,
            D_complex=D_complex,
            L_complex=L_complex,
            linkers={4: tet_linker, 3: tri_linker},
        )
        tet_topo = tet_topo()

        tet_n_metals = 8
        tet_ratios = self._get_ratios(tet_n_metals)
        for rat in tet_ratios:
            print(rat)
            new_name = (
                f"C_{self.cage_set_dict['corner_name']}_"
                f"{self.cage_set_dict['tetratopic']}_"
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
            new_charge = tet_prop['net_charge']*6 + complex_charge

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
                topology_string=tet_topo_name,
                bb_vertices=new_bb_vertices,
                charge=new_charge,
                free_electron_options=new_free_electron_options
            )
            cages_to_build.append(new_cage)
            print('NOT BUILDING ALL RATIOS CURRENTLY!!!!!!')
            break

        # Tritopic homoleptic cages of all symmetries.
        # Get topology as object.
        tri_topo_name, tri_topo = available_topologies(
            string='m4l4spacer'
        )
        tri_topo = tri_topo()

        tri_n_metals = 4
        tri_ratios = self._get_ratios(tri_n_metals)
        for rat in tri_ratios:
            print(rat)
            new_name = (
                f"C_{self.cage_set_dict['corner_name']}_"
                f"{self.cage_set_dict['tritopic']}_"
                f"{tri_topo_name}_"
                f"d{rat[0]}l{rat[1]}_"
                f"SYMM"
            )
            new_bbs = [D_complex, L_complex, tri_linker]
            new_bb_vertices = {
                D_complex: tri_topo.vertices[:rat[0]],
                L_complex: tri_topo.vertices[rat[0]:rat[0]+rat[1]],
                tri_linker: tri_topo.vertices[tri_n_metals:]
            }
            # Merge linker and complex charges.
            complex_charge = rat[0]*int(D_charge)+rat[1]*int(L_charge)
            new_charge = tri_prop['net_charge']*4 + complex_charge

            print(tri_prop['total_unpaired_e'])
            lig_free_e = [
                int(i)*4 for i in tri_prop['total_unpaired_e']
            ]
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
                topology_string=tri_topo_name,
                bb_vertices=new_bb_vertices,
                charge=new_charge,
                free_electron_options=new_free_electron_options
            )
            cages_to_build.append(new_cage)
            print('NOT BUILDING ALL RATIOS CURRENTLY!!!!!!')
            break

        # Prisms of all symmetries.
        # Get topology as object.
        pri_topo_name, pri_topo = available_topologies(
            string='m6l2l3'
        )
        pri_topo = pri_topo()

        pri_n_metals = 6
        pri_ratios = self._get_ratios(pri_n_metals)
        for rat in pri_ratios:
            new_name = (
                f"HeP_{self.cage_set_dict['corner_name']}_"
                f"{self.cage_set_dict['tritopic']}_"
                f"{self.cage_set_dict['tetratopic']}_"
                f"{pri_topo_name}_"
                f"d{rat[0]}l{rat[1]}_"
                f"SYMM"
            )
            new_bbs = [D_complex, L_complex, tri_linker, tet_linker]
            new_bb_vertices = {
                D_complex: pri_topo.vertices[:rat[0]],
                L_complex: pri_topo.vertices[rat[0]:rat[0]+rat[1]],
                tri_linker: pri_topo.vertices[
                    pri_n_metals:pri_n_metals+2
                ],
                tet_linker: pri_topo.vertices[pri_n_metals+2:]
            }
            # Merge linker and complex charges.
            complex_charge = rat[0]*int(D_charge)+rat[1]*int(L_charge)
            new_charge = tri_prop['net_charge']*2
            new_charge += tet_prop['net_charge']*3 + complex_charge

            print(
                tri_prop['total_unpaired_e'],
                tet_prop['total_unpaired_e']
            )
            tri_free_e = [
                int(i)*2 for i in tri_prop['total_unpaired_e']
            ]
            tet_free_e = [
                int(i)*3 for i in tet_prop['total_unpaired_e']
            ]
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
                topology_string=pri_topo_name,
                bb_vertices=new_bb_vertices,
                charge=new_charge,
                free_electron_options=new_free_electron_options
            )
            cages_to_build.append(new_cage)
            print('NOT BUILDING ALL RATIOS CURRENTLY!!!!!!')
            break

        return cages_to_build

    def plot_min_OPs_avg_PV(self, X, Y, T):
        topo_c_m = {
            'm4l4spacer': ('#E074AF', 'o', r'M$_4$L$_4$'),
            'm8l6face': ('#AFE074', 'X', r'M$_8$L$_6$'),
            'm6l2l3': ('#74AFE0', 'P', r'M$_6$L$^a_2$L$^b_3$')
        }

        Cs = [
            topo_c_m[C.topology_string][0]
            for C in self.cages_to_build
        ]
        print(len(X), len(Y), len(Cs))

        fig, ax = plt.subplots(figsize=(8, 5))
        for x, y, t in zip(X, Y, T):
            ax.scatter(
                x,
                y,
                c=topo_c_m[t][0],
                edgecolors='k',
                marker=topo_c_m[t][1],
                alpha=1.0,
                s=80
            )
        # Set number of ticks for x-axis
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel(r'pore volume [$\mathrm{\AA}^3$]', fontsize=16)
        ax.set_ylabel(r'min. $q_{\mathrm{oct}}$', fontsize=16)
        ax.set_xlim(0, 2000)
        ax.set_ylim(0, 1)

        # Implement legend.
        for i in topo_c_m:
            ax.scatter(
                -1000,
                -1000,
                c=topo_c_m[i][0],
                edgecolors='k',
                marker=topo_c_m[i][1],
                alpha=1.0,
                s=80,
                label=topo_c_m[i][2]
            )

        ax.legend(fontsize=16)

        fig.tight_layout()
        fig.savefig(
            f'{self.name}_minOPsVSporevol.pdf',
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()
