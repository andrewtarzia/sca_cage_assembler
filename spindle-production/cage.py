"""Modules defining and building the Cage class and subclasses."""

import json
from os import system
from os.path import exists

import env_set
import pywindow as pw
import stk
import stko
from cage_building import metal_FFs


class UnexpectedNumLigandsError(Exception): ...


class PrecursorNotOptimizedError(Exception): ...


class CageNotOptimizedError(Exception): ...


class Cage:
    """Generic class that builds and analyses stk.ConstructuedMolecules."""

    def __init__(
        self,
        name,
        base_name,
        topology_fn,
        topology_string,
        symmetry_string,
        building_blocks,
        vertex_alignments,
        charge,
        free_electron_options,
        cage_set_dict,
    ):
        self.name = name
        self.base_name = base_name
        self.topology_fn = topology_fn
        self.topology_string = topology_string
        self.symmetry_string = symmetry_string
        self.building_blocks = building_blocks
        self.vertex_alignments = vertex_alignments
        self.topology_graph = self.topology_fn(
            building_blocks=self.building_blocks,
            vertex_alignments=self.vertex_alignments,
        )
        self.unopt_file = f"{self.name}_unopt"
        self.crush_file = f"{self.name}_cru"
        self.uff4mof_file = f"{self.name}_uff"
        self.uff4mof_CG_file = f"{self.name}_uffCG"
        self.uffMD_file = f"{self.name}_prextb"
        self.opt_file = f"{self.name}_optc"
        self.pw_file = f"{self.name}_pw"
        self.op_file = f"{self.name}_OP"
        self.cl_file = f"{self.name}_cl"
        self.sh_file = f"{self.name}_sh"
        self.ls_file = f"{self.name}_LSE"
        self.charge = charge
        self.free_electron_options = free_electron_options
        self.cage_set_dict = cage_set_dict
        self.optimized = None

    def build(self, output_dir):
        cage = stk.ConstructedMolecule(self.topology_graph)
        if not (output_dir / f"{self.unopt_file}.mol").exists():
            cage.write(output_dir / f"{self.unopt_file}.mol")

        self.cage = cage

    def optimize(
        self,
        free_e,
        step_size,
        distance_cut,
        scale_steps,
        output_dir,
    ):
        custom_metal_FFs = metal_FFs(CN=6)

        # Skip if _opt.mol exists.
        if exists(output_dir / f"{self.opt_file}.mol"):
            self.cage = self.cage.with_structure_from_file(
                str(output_dir / f"{self.opt_file}.mol")
            )
            self.optimized = True
            return
        print(f"....optimizing {self.name}")
        self.optimized = None

        # Run if crush output does not exist.
        if not exists(output_dir / f"{self.crush_file}.mol"):
            print(f"..doing collapser optimisation of {self.name}")
            calc_dir = output_dir / f"cage_opt_{self.name}_coll"
            optimizer = stko.Collapser(
                output_dir=calc_dir,
                step_size=step_size,
                distance_cut=distance_cut,
                scale_steps=scale_steps,
            )
            self.cage = optimizer.optimize(mol=self.cage)
            self.cage.write(str(output_dir / f"{self.crush_file}.mol"))
        else:
            self.cage = self.cage.with_structure_from_file(
                str(output_dir / f"{self.crush_file}.mol")
            )

        # Run if uff4mof opt output does not exist.
        if not exists(output_dir / f"{self.uff4mof_CG_file}.mol"):
            cg = True
            maxcyc = 1000
            metal_ligand_bond_order = ""
            calc_dir = (
                output_dir / f"cage_opt_{self.name}_uff"
                if cg is False
                else f"cage_opt_{self.name}_uffCG"
            )
            print(f"..doing UFF4MOF optimisation of {self.name}")
            print(f"Conjugate Gradient: {cg}, Max steps: {maxcyc}")
            gulp_opt = stko.GulpUFFOptimizer(
                gulp_path=env_set.gulp_path(),
                maxcyc=maxcyc,
                metal_FF=custom_metal_FFs,
                metal_ligand_bond_order=metal_ligand_bond_order,
                output_dir=calc_dir,
                conjugate_gradient=cg,
            )
            gulp_opt.assign_FF(self.cage)
            self.cage = gulp_opt.optimize(mol=self.cage)
            self.cage.write(output_dir / f"{self.uff4mof_CG_file}.mol")
        else:
            self.cage = self.cage.with_structure_from_file(
                str(output_dir / f"{self.uff4mof_CG_file}.mol")
            )

        raise SystemExit
        # Run if uff4mof opt output does not exist.
        if not exists(output_dir / f"{self.uff4mof_file}.mol"):
            cg = False
            maxcyc = 1000
            metal_ligand_bond_order = ""
            calc_dir = (
                output_dir / f"cage_opt_{self.name}_uff"
                if cg is False
                else f"cage_opt_{self.name}_uffCG"
            )
            print(f"..doing UFF4MOF optimisation of {self.name}")
            print(f"Conjugate Gradient: {cg}, Max steps: {maxcyc}")
            gulp_opt = stko.GulpUFFOptimizer(
                gulp_path=env_set.gulp_path(),
                maxcyc=maxcyc,
                metal_FF=custom_metal_FFs,
                metal_ligand_bond_order=metal_ligand_bond_order,
                output_dir=calc_dir,
                conjugate_gradient=cg,
            )
            gulp_opt.assign_FF(self.cage)
            self.cage = gulp_opt.optimize(mol=self.cage)
            self.cage.write(output_dir / f"{self.uff4mof_file}.mol")
        else:
            self.cage = self.cage.with_structure_from_file(
                str(output_dir / f"{self.uff4mof_file}.mol")
            )

        raise SystemExit
        # Run if uff4mof MD output does not exist.
        if not exists(output_dir / f"{self.uffMD_file}.mol"):
            print(f"..doing UFF4MOF MD of {self.name}")
            gulp_MD = stko.GulpUFFMDOptimizer(
                gulp_path=env_set.gulp_path(),
                metal_FF=custom_metal_FFs,
                metal_ligand_bond_order="half",
                output_dir=output_dir / f"cage_opt_{self.name}_MD",
                integrator="leapfrog verlet",
                ensemble="nvt",
                temperature=400,
                equilbration=0.1,
                production=2,
                timestep=0.5,
                N_conformers=10,
                opt_conformers=False,
                save_conformers=False,
            )
            gulp_MD.assign_FF(self.cage)
            self.cage = gulp_MD.optimize(self.cage)
            self.cage.write(output_dir / f"{self.uffMD_file}.mol")
        else:
            self.cage = self.cage.with_structure_from_file(
                str(output_dir / f"{self.uffMD_file}.mol")
            )

        raise SystemExit
        try:
            print(f"..........doing XTB optimisation of {self.name}")
            xtb_opt = stko.XTB(
                xtb_path=env_set.xtb_path(),
                output_dir=output_dir / f"cage_opt_{self.name}_xtb",
                gfn_version=2,
                num_cores=6,
                opt_level="normal",
                charge=self.charge,
                num_unpaired_electrons=free_e,
                max_runs=1,
                electronic_temperature=300,
                calculate_hessian=False,
                unlimited_memory=True,
                solvent=self.cage_set_dict["solvent"],
            )
            self.cage = xtb_opt.optimize(mol=self.cage)
            self.cage.write(output_dir / f"{self.opt_file}.mol")
            self.optimized = True
        except (stko.XTBConvergenceError, stko.XTBOptimizerError):
            # Check if the optimisation was even attempted.
            opt_output_file = (
                output_dir
                / f"cage_opt_{self.name}_xtb"
                / "optimization_1.output"
            )
            if not exists(opt_output_file):
                # If not, raise error and exit.
                msg = (
                    "xTB optimisation of cage not even attempted for"
                    f"{self.name}. Try rerunning and check xTB is "
                    "installed."
                )
                raise CageNotOptimizedError(msg)

            # Check if that output file actually contains some
            # steps because xtb may fail on initialisation.
            steps_lines = []
            has_steps = False
            with open(opt_output_file) as f:
                for line in f.readlines():
                    if " CYCLE " in line:
                        steps_lines.append(line)
            if len(steps_lines) > 0:
                has_steps = True
            if has_steps:
                # Set optimized to False, this avoids all analysis.
                self.optimized = False
            else:
                # If not, raise error and exit.
                msg = (
                    "xTB optimisation of cage not even attempted "
                    f"for {self.name}. Try rerunning and check xTB"
                    " is installed."
                )
                raise CageNotOptimizedError(msg)

    def analyze_ligand_strain(
        self,
        metal_atom_no,
        expected_ligands,
        free_e,
    ):
        """Analyse cage ligand geometry for strain."""
        print(f"....analyzing ligand geometry of {self.name}")
        # Collect the atomic positions of the organic linkers in the
        # cage for analysis.
        org_ligs, smiles_keys = get_organic_linkers(
            cage=self.cage,
            metal_atom_nos=(metal_atom_no,),
            file_prefix=f"{self.name}_sg",
        )

        num_unique_ligands = len(set(smiles_keys.values()))
        if num_unique_ligands != expected_ligands:
            msg = (
                f"{self.name} had {num_unique_ligands} unique ligands"
                f", {expected_ligands} were expected. Suggests bad "
                "optimization. Recommend reoptimising structure."
            )
            raise UnexpectedNumLigandsError(msg)

        get_lowest_energy_conformers(
            org_ligs=org_ligs,
            smiles_keys=smiles_keys,
            file_prefix=f"{self.base_name}_sg",
            settings=env_set.crest_conformer_settings(
                solvent=self.cage_set_dict["solvent"],
            ),
        )

        self.calculate_formation_energy(
            org_ligs=org_ligs,
            smiles_keys=smiles_keys,
            file_prefix=f"{self.base_name}_sg",
            cage_free_e=free_e,
        )

        lse_dict = calculate_ligand_SE(
            org_ligs=org_ligs,
            smiles_keys=smiles_keys,
            output_json=f"{self.ls_file}.json",
            file_prefix=f"{self.base_name}_sg",
            solvent=self.cage_set_dict["solvent"],
        )

        imine_torsion_dict = calculate_abs_imine_torsions(
            org_ligs=org_ligs,
        )
        for ol in imine_torsion_dict:
            if len(imine_torsion_dict[ol]) != 4:
                msg = (
                    f"{len(imine_torsion_dict[ol])} minies found, "
                    "but 4 expected."
                )
                raise ValueError(msg)
        planarity_dict = calculate_ligand_planarities(org_ligs=org_ligs)

        self.ls_data = {
            "strain_energies": lse_dict,
            "imine_torsions": imine_torsion_dict,
            "core_planarities": planarity_dict,
        }

    def analyze_metal_strain(self):
        """Analyse cage geometry using order parameters."""
        # Get metal-ligand binder atom bond length.
        print(
            "Warning! Calculate metal-ligand distance is hard-coded for Zn-N"
        )
        self.bl_data = calculate_metal_ligand_distance(
            mol=self.cage,
            metal_atomic_number=30,
            ligand_atomic_number=7,
        )

    def analyze_porosity(self, dump_molecule=False):
        """Analyse cage porosity with pywindow."""
        # Check if output file exists.
        if not exists(f"{self.pw_file}.json"):
            print(f"....analyzing porosity of {self.name}")
            # Load cage into pywindow.
            self.cage.write("temp.xyz")
            pw_cage = pw.MolecularSystem.load_file("temp.xyz")
            pw_cage_mol = pw_cage.system_to_molecule()
            system("rm temp.xyz")

            # Calculate pore size.
            try:
                pw_cage_mol.calculate_pore_diameter_opt()
                pw_cage_mol.calculate_pore_volume_opt()
            except ValueError:
                # Handle failure.
                pw_cage_mol.properties["pore_volume_opt"] = 0
                pw_cage_mol.properties["pore_diameter_opt"] = {
                    "diameter": 0,
                    "atom_1": 0,
                    "centre_of_mass": [0, 0, 0],
                }

            # Save files.
            pw_cage_mol.dump_properties_json(f"{self.pw_file}.json")
            if dump_molecule:
                pw_cage_mol.dump_molecule(
                    f"{self.pw_file}.pdb", include_coms=True
                )

        # Get data.
        with open(f"{self.pw_file}.json") as f:
            self.pw_data = json.load(f)

    def __str__(self):
        return (
            f"{self.__class__.__name__}"
            f"(name={self.name}, topology={self.topology_graph})"
        )

    def __repr__(self):
        return str(self)
