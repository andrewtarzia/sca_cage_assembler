"""Script to build HoCube library."""

import json
import logging
import pathlib
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import stk
import stko
from rdkit.Chem import AllChem as rdkit  # noqa: N813
from utilities import read_lib

logger = logging.getLogger(__name__)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)


def unit_vector(vector: np.ndarray) -> np.ndarray:
    """Returns the unit vector of the vector."""
    return vector / np.linalg.norm(vector)


def vector_angle(vector1: np.ndarray, vector2: np.ndarray) -> float:
    """Returns the angle between two vectors in radians.

    From: https://stackoverflow.com/a/13849249

    Parameters:
        vector1:
            The first vector.

        vector2:
            The second vector.

    Returns:
        The angle between `vector1` and `vector2` in radians.

    """
    v1_u = unit_vector(vector1)
    v2_u = unit_vector(vector2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


class GeometryAnalyser:
    """Tools for analysing the geometry of molecules.

    .. warning::
        This code is only present in the latest versions of stko
        that require Python 3.11!

    """

    def _get_paths(
        self,
        molecule: stk.Molecule,
        path_length: int,
    ) -> tuple[tuple[int, ...], ...]:
        return rdkit.FindAllPathsOfLengthN(
            mol=molecule.to_rdkit_mol(),
            length=path_length,
            useBonds=False,
            useHs=True,
        )

    def calculate_bonds(
        self,
        molecule: stk.Molecule,
    ) -> dict[tuple[str, str], list[float]]:
        """Calculate bond lengths for all `stk.Molecule.get_bonds()`.

        Parameters:
            molecule:
                The molecule to analyse.

        Returns:
            Dictionary of bonds organised by element pair.

        """
        position_matrix = molecule.get_position_matrix()
        lengths: dict[tuple[str, str], list[float]] = defaultdict(list)
        for bond in molecule.get_bonds():
            a1id = bond.get_atom1().get_id()
            a2id = bond.get_atom2().get_id()
            a, b = sorted(
                (
                    bond.get_atom1().__class__.__name__,
                    bond.get_atom2().__class__.__name__,
                )
            )
            lengths[(a, b)].append(
                stko.get_atom_distance(position_matrix, a1id, a2id)
            )

        return lengths

    def calculate_angles(
        self,
        molecule: stk.Molecule,
    ) -> dict[tuple[str, str, str], list[float]]:
        """Calculate angles for all angles defined by molecule bonding.

        Parameters:
            molecule:
                The molecule to analyse.

        Returns:
            Dictionary of angles organised by element triplet.

        """
        position_matrix = molecule.get_position_matrix()
        angles: dict[tuple[str, str, str], list[float]] = defaultdict(list)
        for a_ids in self._get_paths(molecule, 3):
            atoms = list(molecule.get_atoms(atom_ids=a_ids))
            atom1 = atoms[0]
            atom2 = atoms[1]
            atom3 = atoms[2]
            angle_type_option1 = (
                atom1.__class__.__name__,
                atom2.__class__.__name__,
                atom3.__class__.__name__,
            )
            angle_type_option2 = (
                atom3.__class__.__name__,
                atom2.__class__.__name__,
                atom1.__class__.__name__,
            )

            vector1 = (
                position_matrix[atom2.get_id()]
                - position_matrix[atom1.get_id()]
            )
            vector2 = (
                position_matrix[atom2.get_id()]
                - position_matrix[atom3.get_id()]
            )

            if angle_type_option1 in angles:
                angles[angle_type_option1].append(
                    np.degrees(vector_angle(vector1, vector2))
                )
            elif angle_type_option2 in angles:
                angles[angle_type_option2].append(
                    np.degrees(vector_angle(vector1, vector2))
                )
            else:
                angles[angle_type_option1].append(
                    np.degrees(vector_angle(vector1, vector2))
                )

        return angles

    def calculate_torsions(
        self,
        molecule: stk.Molecule,
    ) -> dict[tuple[str, ...], list[float]]:
        """Calculate all torsions defined by molecule bonding.

        Parameters:
            molecule:
                The molecule to analyse.

        Returns:
            Dictionary of torsions organised by elements.

        """
        position_matrix = molecule.get_position_matrix()

        torsions: dict[tuple[str, ...], list[float]] = defaultdict(list)
        for a_ids in self._get_paths(molecule, 4):
            atoms = list(molecule.get_atoms(atom_ids=a_ids))
            atom1 = atoms[0]
            atom2 = atoms[1]
            atom3 = atoms[2]
            atom4 = atoms[3]
            torsion_type_option1 = (
                atom1.__class__.__name__,
                atom2.__class__.__name__,
                atom3.__class__.__name__,
                atom4.__class__.__name__,
            )
            torsion_type_option2 = (
                atom4.__class__.__name__,
                atom3.__class__.__name__,
                atom2.__class__.__name__,
                atom1.__class__.__name__,
            )

            if torsion_type_option1 in torsions:
                torsions[torsion_type_option1].append(
                    stko.calculate_dihedral(
                        pt1=position_matrix[atom1.get_id()],
                        pt2=position_matrix[atom2.get_id()],
                        pt3=position_matrix[atom3.get_id()],
                        pt4=position_matrix[atom4.get_id()],
                    )
                )
            elif torsion_type_option2 in torsions:
                torsions[torsion_type_option2].append(
                    stko.calculate_dihedral(
                        pt1=position_matrix[atom4.get_id()],
                        pt2=position_matrix[atom3.get_id()],
                        pt3=position_matrix[atom2.get_id()],
                        pt4=position_matrix[atom1.get_id()],
                    )
                )
            else:
                torsions[torsion_type_option1].append(
                    stko.calculate_dihedral(
                        pt1=position_matrix[atom1.get_id()],
                        pt2=position_matrix[atom2.get_id()],
                        pt3=position_matrix[atom3.get_id()],
                        pt4=position_matrix[atom4.get_id()],
                    )
                )

        return torsions


def parity_energies(
    cage_set_lib,
    cage_directory,
    old_cage_directory,
    figure_directory,
) -> None:
    """Plot energy parity."""
    fig, (ax, ax1) = plt.subplots(ncols=2, figsize=(10, 5))

    name_convention = {
        "cl1_quad2_8": "quad2_8",
        "cl1_quad2_12": "quad2_12",
        "cl1_quad2_2": "quad2_2",
        "cl1_quad2_3": "quad2_3",
        "cl1_quad2_9": "quad2_9",
        "cl1_quad2_10": "quad2_10",
        "cl1_quad2_5": "quad2_5",
        "cl1_quad2_16": "quad2_16",
        "cl1_quad2_17": "quad2_17",
    }

    for name in cage_set_lib:
        ey_files = sorted(cage_directory.glob(f"*{name}*_optc.ey"))
        energies = {}
        for eyf in ey_files:
            with eyf.open("r") as f:
                data = f.readlines()[0]
            energies[eyf.stem] = float(data)

        ey_files = sorted(old_cage_directory.glob(f"*{name}*_optc.ey"))
        old_energies = {}
        for eyf in ey_files:
            with eyf.open("r") as f:
                data = f.readlines()[0]
            old_energies[eyf.stem] = float(data)

        paired_keys = [i for i in energies if i in old_energies]
        if len(paired_keys) == 0:
            continue
        min_energy = min(energies.values())
        min_old_energy = min(old_energies.values())
        ax.scatter(
            [(old_energies[i] - min_old_energy) * 2625.5 for i in paired_keys],
            [(energies[i] - min_energy) * 2625.5 for i in paired_keys],
            edgecolors="k",
            marker="o",
            alpha=1.0,
            s=40,
            label=name_convention[name],
        )
        ax1.scatter(
            [(old_energies[i] - min_old_energy) * 2625.5 for i in paired_keys],
            [(energies[i] - min_energy) * 2625.5 for i in paired_keys],
            edgecolors="k",
            marker="o",
            alpha=1.0,
            s=80,
        )

    # Set number of ticks for x-axis
    ax.tick_params(axis="both", which="major", labelsize=16)
    ax.set_xlabel("2022 rel. GFN2-xTB energy [kJmol$^{-1}$]", fontsize=16)
    ax.set_ylabel("2025 rel. GFN2-xTB energy [kJmol$^{-1}$]", fontsize=16)
    ax.set_xlim(0, 1000)
    ax.set_ylim(0, 1000)
    ax.plot((0, 1000), (0, 1000), c="k", zorder=-2)
    ax.legend(fontsize=16)

    ax1.tick_params(axis="both", which="major", labelsize=16)
    ax1.set_xlabel("2022 rel. GFN2-xTB energy [kJmol$^{-1}$]", fontsize=16)
    ax1.set_xlim(0, 200)
    ax1.set_ylim(0, 200)
    ax1.plot((0, 200), (0, 200), c="k", zorder=-2)

    fig.tight_layout()
    fig.savefig(
        figure_directory / "parities_energy.pdf",
        dpi=720,
        bbox_inches="tight",
    )
    fig.savefig(
        figure_directory / "parities_energy.png",
        dpi=720,
        bbox_inches="tight",
    )
    plt.close()


def parity_strain_energies(
    cage_set_lib,
    old_cage_directory,
    cage_directory,
    figure_directory,
) -> None:
    """Plot energy parity."""
    fig, (ax, ax1) = plt.subplots(ncols=2, figsize=(10, 5))

    name_convention = {
        "cl1_quad2_8": "quad2_8",
        "cl1_quad2_12": "quad2_12",
        "cl1_quad2_2": "quad2_2",
        "cl1_quad2_3": "quad2_3",
        "cl1_quad2_9": "quad2_9",
        "cl1_quad2_10": "quad2_10",
        "cl1_quad2_5": "quad2_5",
        "cl1_quad2_16": "quad2_16",
        "cl1_quad2_17": "quad2_17",
    }

    for name in cage_set_lib:
        lses = {}
        cs_file = cage_directory / f"{name}_CS.json"
        with cs_file.open("r") as f:
            cs_data = json.load(f)
        for struct, sdata in cs_data.items():
            if not sdata["optimized"]:
                continue

            lses[struct] = sum(
                [
                    sdata["li_prop"]["strain_energies"][i]
                    for i in sdata["li_prop"]["strain_energies"]
                ]
            )

        old_lses = {}
        old_cs_file = old_cage_directory / f"{name}_CS.json"
        try:
            with old_cs_file.open("r") as f:
                old_cs_data = json.load(f)
        except FileNotFoundError:
            continue
        for struct, sdata in old_cs_data.items():
            if not sdata["optimized"]:
                continue

            old_lses[struct] = sum(
                [
                    sdata["li_prop"]["strain_energies"][i]
                    for i in sdata["li_prop"]["strain_energies"]
                ]
            )

        paired_keys = [i for i in lses if i in old_lses]
        if len(paired_keys) == 0:
            continue

        min_energy = min(lses.values())
        min_old_energy = min(old_lses.values())
        ax.scatter(
            [(old_lses[i] - min_old_energy) for i in paired_keys],
            [(lses[i] - min_energy) for i in paired_keys],
            edgecolors="k",
            marker="o",
            alpha=1.0,
            s=40,
            label=name_convention[name],
        )
        ax1.scatter(
            [(old_lses[i] - min_old_energy) for i in paired_keys],
            [(lses[i] - min_energy) for i in paired_keys],
            edgecolors="k",
            marker="o",
            alpha=1.0,
            s=80,
        )

    # Set number of ticks for x-axis
    ax.tick_params(axis="both", which="major", labelsize=16)
    ax.set_xlabel("2022 sum strain energy [kJmol$^{-1}$]", fontsize=16)
    ax.set_ylabel("2025 sum strain energy [kJmol$^{-1}$]", fontsize=16)
    ax.set_xlim(0, 1000)
    ax.set_ylim(0, 1000)
    ax.plot((0, 1000), (0, 1000), c="k", zorder=-2)
    ax.legend(fontsize=16)

    ax1.tick_params(axis="both", which="major", labelsize=16)
    ax1.set_xlabel("2022 sum strain energy [kJmol$^{-1}$]", fontsize=16)
    ax1.set_xlim(0, 200)
    ax1.set_ylim(0, 200)
    ax1.plot((0, 200), (0, 200), c="k", zorder=-2)

    fig.tight_layout()
    fig.savefig(
        figure_directory / "parities_rellse.pdf",
        dpi=720,
        bbox_inches="tight",
    )
    fig.savefig(
        figure_directory / "parities_rellse.png",
        dpi=720,
        bbox_inches="tight",
    )
    plt.close()


def parity_geometries(
    cage_set_lib,
    old_cage_directory,
    cage_directory,
    figure_directory,
) -> None:
    """Plot energy parity."""
    fig, ax = plt.subplots(figsize=(8, 5))

    analyser = GeometryAnalyser()
    n_zn_distances = []
    old_n_zn_distances = []
    for name in cage_set_lib:
        lses = {}

        structures = sorted(cage_directory.glob(f"*{name}*_optc.mol"))
        # geoms = {}
        for sf in structures:
            mol = stk.BuildingBlock.init_from_file(str(sf))
            n_zn_distances.extend(analyser.calculate_bonds(mol)[("N", "Zn")])
            # geoms[sf.stem] = {}
            # geoms[sf.stem]["bonds"] = analyser.calculate_bonds(mol)
            # geoms[sf.stem]["angles"] = analyser.calculate_angles(mol)
            # geoms[sf.stem]["torsions"] = analyser.calculate_torsions(mol)

        structures = sorted(old_cage_directory.glob(f"*{name}*_optc.mol"))
        # old_geoms = {}
        for sf in structures:
            mol = stk.BuildingBlock.init_from_file(str(sf))
            old_n_zn_distances.extend(
                analyser.calculate_bonds(mol)[("N", "Zn")]
            )
            # old_geoms[sf.stem] = {}
            # old_geoms[sf.stem]["bonds"] = analyser.calculate_bonds(mol)
            # old_geoms[sf.stem]["angles"] = analyser.calculate_angles(mol)
            # old_geoms[sf.stem]["torsions"] = analyser.calculate_torsions(mol)

    xmin = 1.8
    xmax = 2.5
    xwidth = 0.01
    xbins = np.arange(xmin - xwidth, xmax + xwidth, xwidth)

    ax.hist(
        x=old_n_zn_distances,
        bins=xbins,
        density=True,
        bottom=0,
        histtype="stepfilled",
        stacked=True,
        linewidth=1.0,
        edgecolor="k",
        label="2022",
    )
    ax.hist(
        x=n_zn_distances,
        bins=xbins,
        density=True,
        bottom=15,
        histtype="stepfilled",
        stacked=True,
        linewidth=1.0,
        edgecolor="k",
        label="2025",
    )

    # Set number of ticks for x-axis
    ax.tick_params(axis="both", which="major", labelsize=16)
    ax.set_xlabel("N-Zn distance [AA]", fontsize=16)
    ax.set_ylabel("density", fontsize=16)
    ax.set_yticks([])
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(
        figure_directory / "parities_geometries.pdf",
        dpi=720,
        bbox_inches="tight",
    )
    fig.savefig(
        figure_directory / "parities_geometries.png",
        dpi=720,
        bbox_inches="tight",
    )
    plt.close()


def main() -> None:
    """Run script."""
    script_directory = pathlib.Path(__file__).parent.resolve()
    data_directory = script_directory / ".." / "data"
    working_dir = pathlib.Path(
        "/home/atarzia/workingspace/spindle_project/rerun_production/"
    )

    cage_set_lib_file = data_directory / "spindle_library.json"
    cage_directory = working_dir / "cage_library"

    figure_directory = working_dir / "figures"
    figure_directory.mkdir(parents=True, exist_ok=True)

    cage_set_lib = read_lib(cage_set_lib_file)

    parity_energies(
        cage_set_lib=cage_set_lib,
        cage_directory=cage_directory,
        old_cage_directory=working_dir / ".." / "ey_files",
        figure_directory=figure_directory,
    )
    parity_strain_energies(
        cage_set_lib=cage_set_lib,
        cage_directory=cage_directory,
        old_cage_directory=working_dir / ".." / "cs_jsons",
        figure_directory=figure_directory,
    )
    parity_geometries(
        cage_set_lib=cage_set_lib,
        cage_directory=cage_directory,
        old_cage_directory=working_dir / ".." / "optc_from_zip_and_clean",
        figure_directory=figure_directory,
    )


if __name__ == "__main__":
    main()
