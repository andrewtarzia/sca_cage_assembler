"""Modules defining and building the Cage class and subclasses."""

import json
from itertools import product
from os.path import join

import stk
from cage import Cage
from cage_building import available_topologies
from facebuildingblock import FaceBuildingBlock
from symmetries import M8L6_Symmetry, M8L6Knot_Symmetry
from utilities import reorient_linker


class CageSet:
    """Class that builds and analyses stk.ConstructuedMolecules."""

    def __init__(
        self,
        name,
        cage_set_dict,
        complex_dicts,
        ligand_dicts,
        ligand_dir,
        complex_dir,
        cage_dir,
    ):
        self.name = name
        self.properties_file = cage_dir / f"{self.name}_CS.json"
        self.measures_file = cage_dir / f"{self.name}_measures.json"
        self.cage_set_dict = cage_set_dict
        self.complex_dicts = complex_dicts
        self.ligand_dicts = ligand_dicts
        self.cages_to_build = self.define_cages_to_build(
            ligand_dir=ligand_dir,
            complex_dir=complex_dir,
        )
        self.built_cage_properties = {}

    def get_sum_lig_strain_energy(self, cage_name):
        C_data = self.built_cage_properties[cage_name]
        if C_data["optimized"] is False:
            return None

        return sum(
            [
                C_data["li_prop"]["strain_energies"][i]
                for i in C_data["li_prop"]["strain_energies"]
            ]
        )

    def get_min_imine_torision(self, cage_name):
        C_data = self.built_cage_properties[cage_name]
        if C_data["optimized"] is False:
            return None

        return min(
            [
                j
                for i in C_data["li_prop"]["imine_torsions"]
                for j in C_data["li_prop"]["imine_torsions"][i]
            ]
        )

    def get_pore_diameter(self, cage_name):
        C_data = self.built_cage_properties[cage_name]
        if C_data["optimized"] is False:
            return None

        return C_data["pw_prop"]["pore_diameter_opt"]["diameter"]

    def get_max_ML_distance(self, cage_name):
        C_data = self.built_cage_properties[cage_name]
        if C_data["optimized"] is False:
            return None

        return max([i for i in C_data["bl_prop"]])

    def define_cages_to_build(self):
        """Defines the name and objects of all cages to build."""
        msg = f"Not implemented for {self.__class__}"
        raise NotImplementedError(msg)

    def load_properties(self):
        """Load class from JSON file."""
        if self.properties_file.exists():
            with self.properties_file.open() as f:
                self.built_cage_properties = json.load(f)
        else:
            msg = f"{self.properties_file} does not exist"
            raise FileNotFoundError(msg)

    def dump_properties(self):
        """Dump class to JSON file."""
        with self.properties_file.open("w") as f:
            json.dump(self.built_cage_properties, f, indent=4)

    def __str__(self):
        return (
            f"{self.__class__.__name__}"
            f"(name={self.name})\n"
            f"{self.cage_set_dict}"
        )

    def __repr__(self):
        return str(self)

    def _get_no_vertices(self, string):
        """Get the number of vertices for a given topology."""
        topologies = {
            "m4l4spacer": 8,
            "m8l6face": 14,
            "m6l2l3": 11,
        }

        try:
            return topologies[string]
        except KeyError:
            raise KeyError(f"{string} not in {topologies.keys()}")

    def _get_complex_info(self, complex_dir):
        D_complex_name = [i for i in self.complex_dicts if "del" in i][0]
        L_complex_name = [i for i in self.complex_dicts if "lam" in i][0]
        L_complex = self._load_complex(
            complex_name=L_complex_name, complex_dir=complex_dir
        )
        D_complex = self._load_complex(
            complex_name=D_complex_name, complex_dir=complex_dir
        )

        return D_complex_name, D_complex, L_complex_name, L_complex

    def _get_complex_properties(self, complex_name):
        charge = self.complex_dicts[complex_name]["total_charge"]
        free_e = self.complex_dicts[complex_name]["unpaired_e"]

        return charge, free_e

    def _get_ligand(self, type_name, ligand_dir):
        prop = self.ligand_dicts[self.cage_set_dict[type_name]]
        linker = self._load_ligand(
            ligand_name=self.cage_set_dict[type_name], ligand_dir=ligand_dir
        )

        return prop, linker

    def _get_rot_vertices(self, string):
        """Get the list of rotatable vertices for a given topology.

        Only ligand vertices are rotatable in this case.

        # TODO: Currently only defined for cube (90 deg). Add tri-face.

        """
        if string in ["m4l4spacer", "m6l2l3"]:
            msg = "Currently only defined for cube (90 deg). Add tri."
            raise NotImplementedError(msg)

        topologies = {
            "m4l4spacer": [4, 5, 6, 7],
            "m8l6face": [8, 9, 10, 11, 12, 13],
            "m6l2l3": [8, 9, 10],
        }

        try:
            return topologies[string]
        except KeyError as e:
            msg = f"{string} not in {topologies.keys()}"
            raise KeyError(msg) from e

    def _get_ratios(self, n_metals):
        rng = range(n_metals + 1)
        rats = []
        for i in product(rng, rng):
            if i[0] + i[1] == n_metals:
                rats.append(i)
        return rats

    def _load_complex(self, complex_name, complex_dir):
        complex = stk.BuildingBlock.init_from_file(
            join(complex_dir, f"{complex_name}_opt.mol"),
            functional_groups=[stk.BromoFactory()],
        )

        return complex

    def _load_ligand(self, ligand_name, ligand_dir):
        ligand = stk.BuildingBlock.init_from_file(
            join(ligand_dir, f"{ligand_name}_opt.mol"),
            functional_groups=[stk.BromoFactory()],
        )

        return ligand

    def get_cage_symmetries(
        self,
        string,
        D_complex,
        L_complex,
        linkers,
    ):
        """Returns cage symmetries for a given topology."""
        if string == "m8l6face":
            symm_list = {}
            linker = linkers[4]

            # Predefined list of symmetries.
            symm_c = M8L6_Symmetry(
                D_complex=D_complex,
                L_complex=L_complex,
                linker=linker,
            )
            symm_list["d2"] = symm_c.d2()
            symm_list["th1"] = symm_c.th1()
            symm_list["th2"] = symm_c.th2()
            symm_list["td"] = symm_c.td()
            symm_list["tl"] = symm_c.tl()
            # symm_list["s41"] = symm_c.s41()
            # symm_list["s42"] = symm_c.s42()
            # symm_list["s61"] = symm_c.s61()
            symm_list["s62"] = symm_c.s62()
            # symm_list["d31"] = symm_c.d31()
            symm_list["d32"] = symm_c.d32()
            symm_list["d31n"] = symm_c.d31n()
            symm_list["d32n"] = symm_c.d32n()
            # symm_list["c2v"] = symm_c.c2v()
            # symm_list["c2h"] = symm_c.c2h()

        elif string == "m8l6knot":
            symm_list = {}
            linker = linkers[4]

            # Predefined list of symmetries.
            symm_c = M8L6Knot_Symmetry(
                D_complex=D_complex,
                L_complex=L_complex,
                linker=linker,
            )
            symm_list["d3c3"] = symm_c.d3c3()
        else:
            msg = f"{string} not in defined"
            raise KeyError(msg)

        return symm_list

    def iterate_over_symmetries(
        self,
        base_name,
        topo_name,
        topo_fn,
        symmetries_to_build,
        charge_prop,
        mult_prop,
    ):
        """Iterates over symmetry options and defines .Cage."""
        cages_to_build = []

        for name_string in symmetries_to_build:
            new_name = f"{base_name}_{name_string}"
            building_blocks = symmetries_to_build[name_string][
                "building_blocks"
            ]
            vertex_alignments = symmetries_to_build[name_string][
                "vertex_alignments"
            ]
            rat = symmetries_to_build[name_string]["ratio"]

            # Merge linker and complex charges.
            complex_charge = rat[0] * charge_prop["D"]
            complex_charge += rat[1] * charge_prop["L"]
            new_charge = charge_prop["4"] + charge_prop["3"] + complex_charge

            compl_free_e = [
                int(i) * rat[0] + int(j) * rat[1]
                for i, j in zip(mult_prop["D"], mult_prop["L"])
            ]
            new_free_electron_options = []
            for opt in product(
                mult_prop["3"],
                mult_prop["4"],
                compl_free_e,
            ):
                new_free_electron_options.append(opt[0] + opt[1] + opt[2])

            new_cage = Cage(
                name=new_name,
                base_name=base_name,
                topology_fn=topo_fn,
                building_blocks=building_blocks,
                vertex_alignments=vertex_alignments,
                topology_string=topo_name,
                symmetry_string=name_string,
                charge=new_charge,
                free_electron_options=new_free_electron_options,
                cage_set_dict=self.cage_set_dict,
            )
            cages_to_build.append(new_cage)

        return cages_to_build


class HoCube(CageSet):
    """Class that builds and analyses stk.ConstructuedMolecules.

    Represents homoleptic cube cages with all necessary symmetries
    and orientations.

    """

    def _load_ligand(self, ligand_name, ligand_dir):
        return FaceBuildingBlock.init_from_file(
            str(ligand_dir / f"{ligand_name}_opt.mol"),
            functional_groups=[stk.BromoFactory()],
        )

    def _get_ligand(self, type_name, ligand_dir):
        prop = self.ligand_dicts[self.cage_set_dict[type_name]]
        linker = self._load_ligand(
            ligand_name=self.cage_set_dict[type_name],
            ligand_dir=ligand_dir,
        )

        temp_linker = reorient_linker(linker)

        # Set functional group ordering based on long axis.
        fg_centroids = tuple(
            temp_linker.get_centroid(
                atom_ids=fg.get_placer_ids(),
            )
            for fg in temp_linker.get_functional_groups()
        )
        plus_minus_fg_id = next(
            i
            for i, cent in enumerate(fg_centroids)
            if cent[0] > 0 and cent[1] < 0
        )
        fg1_id = plus_minus_fg_id
        fg2_id, fg3_id, fg4_id = tuple(
            i
            for i in range(temp_linker.get_num_functional_groups())
            if i != fg1_id
        )
        new_fgs = tuple(temp_linker.get_functional_groups())
        linker = temp_linker.with_functional_groups(
            functional_groups=(
                new_fgs[fg1_id],
                new_fgs[fg2_id],
                new_fgs[fg3_id],
                new_fgs[fg4_id],
            )
        )
        return prop, linker

    def define_cages_to_build(self, ligand_dir, complex_dir):
        """Defines the name and objects of all cages to build."""
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
            type_name="tetratopic",
            ligand_dir=ligand_dir,
        )

        # Get topology function as object to be used in following list.
        # Homoleptic cage with tetratopic ligand.
        tet_topo_names = ["m8l6face", "m8l6knot"]
        cages_to_build = []
        for topo_name in tet_topo_names:
            tet_topo_fn = available_topologies(string=topo_name)

            symmetries_to_build = self.get_cage_symmetries(
                string=topo_name,
                D_complex=D_complex,
                L_complex=L_complex,
                linkers={4: tet_linker},
            )

            temp_cages_to_build = self.iterate_over_symmetries(
                base_name=(
                    f"C_{self.cage_set_dict['corner_name']}_"
                    f"{self.cage_set_dict['tetratopic']}"
                ),
                topo_name=topo_name,
                topo_fn=tet_topo_fn,
                symmetries_to_build=symmetries_to_build,
                # Set charge properties based on ligand occurances.
                charge_prop={
                    "D": int(D_charge),
                    "L": int(L_charge),
                    "3": 0,
                    "4": tet_prop["net_charge"] * 6,
                },
                # Set free e properties based on ligand occurances.
                mult_prop={
                    "D": D_free_e,
                    "L": L_free_e,
                    "3": [0],
                    "4": [int(i) * 6 for i in tet_prop["total_unpaired_e"]],
                },
            )

            for i in temp_cages_to_build:
                cages_to_build.append(i)  # noqa: PERF402

        return cages_to_build
