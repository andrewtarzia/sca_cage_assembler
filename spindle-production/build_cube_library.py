"""Script to build HoCube library."""

import logging
import pathlib

import stko
from cage_set import HoCube
from utilities import read_lib

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)


def build_cages(
    ligands,
    complexes,
    cage_set_lib,
    ligand_directory,
    complex_directory,
    cage_directory,
):
    cage_sets = []
    for name in cage_set_lib:
        cage_set_c = cage_set_lib[name]
        compl_names = cage_set_c["corners"]
        comps = {i: complexes[i] for i in compl_names}

        cage_set = HoCube(
            name=name,
            cage_set_dict=cage_set_c,
            complex_dicts=comps,
            ligand_dicts=ligands,
            ligand_dir=ligand_directory,
            complex_dir=complex_directory,
            cage_dir=cage_directory,
        )

        if cage_set.properties_file.exists():
            cage_set.load_properties()

        for cage in cage_set.cages_to_build:
            logging.info("building cage %s...", cage.name)
            cage.build(output_dir=cage_directory)

            default_free_e = cage.free_electron_options[0]
            # Use a slightly different collapser threshold for
            # different topologies.
            if cage.topology_string in ["m8l6face"]:
                step_size = 0.05
                distance_cut = 2.5
                scale_steps = False
                expected_ligands = 1
                coll_fun = stko.Collapser
            elif cage.topology_string in ["m8l6knot"]:
                step_size = 0.1
                distance_cut = 2
                scale_steps = True
                expected_ligands = 1
                coll_fun = stko.CollapserMC
            else:
                raise NotImplementedError

            # Check if structure has previously had optimisation
            # attempted with failure.
            try:
                if (
                    cage_set.built_cage_properties[cage.name]["optimized"]
                    is False
                ):
                    cage.optimized = False

            except KeyError:
                pass

            if cage.optimized is None:
                # Run optimisation - which handles successfully
                # completed molecules.
                cage.optimize(
                    free_e=default_free_e,
                    step_size=step_size,
                    distance_cut=distance_cut,
                    scale_steps=scale_steps,
                    output_dir=cage_directory,
                    coll_fun=coll_fun,
                )

            if cage.optimized is False:
                cage_set.built_cage_properties[cage.name] = {
                    "optimized": cage.optimized,
                }

            else:
                cage.analyze_ligand_strain(
                    # Assumes only one type of metal atom.
                    metal_atom_no=next(
                        cage_set.complex_dicts[i]["metal_atom_no"]
                        for i in cage_set.complex_dicts
                    ),
                    expected_ligands=expected_ligands,
                    free_e=default_free_e,
                    output_dir=cage_directory,
                )
                cage.analyze_energy(
                    output_dir=cage_directory,
                    free_e=default_free_e,
                )

                cage.analyze_metal_strain()
                cage.analyze_porosity(output_dir=cage_directory)
                cage_set.built_cage_properties[cage.name] = {
                    "optimized": cage.optimized,
                    "pw_prop": cage.pw_data,
                    "li_prop": cage.ls_data,
                    "bl_prop": cage.bl_data,
                }

            # Dump to JSON.
            cage_set.dump_properties()

        cage_sets.append(cage_set)

    return cage_sets


def main() -> None:
    """Run script."""
    script_directory = pathlib.Path(__file__).parent.resolve()
    data_directory = script_directory / ".." / "data"
    working_dir = pathlib.Path(
        "/home/atarzia/onbear/tarziaa-cont1/local/spindle_project/rerun_production/"
    )

    ligand_lib_file = data_directory / "spindle_ligand_library.json"
    complex_lib_file = data_directory / "cube_complex_library.json"
    cage_set_lib_file = data_directory / "spindle_library.json"
    ligand_directory = working_dir / "ligand_library"
    ligand_directory.mkdir(exist_ok=True, parents=True)
    complex_directory = working_dir / "complex_library"
    complex_directory.mkdir(exist_ok=True, parents=True)
    cage_directory = working_dir / "cage_library"
    cage_directory.mkdir(exist_ok=True, parents=True)

    cage_set_lib = read_lib(cage_set_lib_file)
    complexes = read_lib(complex_lib_file)
    ligands = read_lib(ligand_lib_file)

    # Build and optimise all organic molecules in lib.
    build_cages(
        ligands=ligands,
        complexes=complexes,
        cage_set_lib=cage_set_lib,
        ligand_directory=ligand_directory,
        complex_directory=complex_directory,
        cage_directory=cage_directory,
    )


if __name__ == "__main__":
    main()
