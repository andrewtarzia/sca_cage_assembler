#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Utilitiy functions for setting up CP2K DFT calculations.

Author: Andrew Tarzia

Date Created: 20 Jan 2022

"""

import os

import stk


class CP2KCalculator:
    """

    """

    def __init__(self, job_name):
        self._job_name = job_name

    def write_slurm_file(self, output_directory, hours):

        string = (
            '#!/bin/bash\n\n'
            f'#SBATCH --job-name={self._job_name}\n'
            '#SBATCH --nodes=10\n'
            '#SBATCH --tasks-per-node=128\n'
            '#SBATCH --cpus-per-task=1\n'
            f'#SBATCH --time={hours}:00:00\n\n'
            '#SBATCH --account=e05-discov-kim\n'
            '#SBATCH --partition=standard\n'
            '#SBATCH --qos=standard\n\n'
            '# Load the relevent CP2K module\n'
            'module load cp2k\n'
            'export OMP_NUM_THREADS=1\n\n'
            'srun --hint=nomultithread --distribution=block:block '
            f'cp2k.popt -i {self._job_name}.in -o {self._job_name}.out'
        )

        slurm_output = os.path.join(
            output_directory, f'{self._job_name}.slurm'
        )
        with open(slurm_output, 'w') as f:
            f.write(string)

    def _global_section(self):
        raise NotImplementedError()

    def _force_section(
        self,
        charge,
        cutoff,
        rel_cutoff,
        guess=None,
        solvent=None,
    ):

        if solvent is None:
            solvent_section = ''
        else:
            solvent_section = (
                '  &SCCS\n'
                f'    RELATIVE_PERMITTIVITY {solvent}\n'
                '  &END SCCS\n'
            )

        if guess is None:
            guess = 'ATOMIC'

        string = (
'! Force Evaluation\n'
'! ----------------\n'
'&FORCE_EVAL\n'
'  METHOD QS\n'
'  &DFT\n'
'    BASIS_SET_FILE_NAME BASIS_MOLOPT\n'
'    POTENTIAL_FILE_NAME GTH_POTENTIALS\n'
f'    CHARGE {charge}\n'
'  &QS\n'
'   EPS_DEFAULT 1.0E-10\n'
'  &END QS\n'
f'{solvent_section}'
'  &MGRID\n'
'    NGRIDS 5\n'
f'   CUTOFF {cutoff}\n'
f'   REL_CUTOFF {rel_cutoff}\n'
'  &END MGRID\n'
'  &SCF\n'
f'      SCF_GUESS {guess}\n'
'      EPS_SCF 5.0E-7\n'
'      MAX_SCF 300\n'
'      &OT\n'
'        PRECONDITIONER FULL_SINGLE_INVERSE\n'
'        MINIMIZER DIIS\n'
'      &END OT\n'
'      &OUTER_SCF ! Repeat the inner SCF cycle 10 times\n'
'        MAX_SCF 10\n'
'        EPS_SCF 5.0E-7 ! Must match the above\n'
'      &END\n'
'      &PRINT\n'
'        &RESTART ON\n'
f'        FILENAME {self._job_name}.res\n'
'        &END\n'
'      &END\n'
'  &END SCF\n'
'  ! Exchange-Correlation\n'
'  ! -------------------\n'
'  &XC\n'
'  &VDW_POTENTIAL\n'
'    POTENTIAL_TYPE PAIR_POTENTIAL\n'
'    &PAIR_POTENTIAL\n'
'      TYPE DFTD3\n'
'      CALCULATE_C9_TERM .TRUE.\n'
'      PARAMETER_FILE_NAME dftd3.dat\n'
'      REFERENCE_FUNCTIONAL PBE\n'
'      R_CUTOFF 10.\n'
'      EPS_CN 0.01\n'
'    &END PAIR_POTENTIAL\n'
'  &END VDW_POTENTIAL\n'
'  &XC_FUNCTIONAL PBE\n'
'  &END XC_FUNCTIONAL\n'
'&END XC\n'
'&END DFT\n'
'  &PRINT\n'
'    &FORCES\n'
'    &END\n'
'  &END\n'
        )

        return string

    def _cell_section(self, molecule):

        val = round(molecule.get_maximum_diameter()+10, 2)

        string = (
            '&SUBSYS\n\n'
            '! Cell\n'
            '! ----\n'
            '&CELL\n'
            f'  ABC {val} {val} {val}\n'
            '  ALPHA_BETA_GAMMA 90.0 90.0 90.0\n'
            '  PERIODIC NONE\n'
            '&END CELL\n\n'
        )
        return string

    def _coordinates_section(self, molecule):
        string = (
            '! Coordinates\n'
            '! -----------\n'
            '&COORD\n'
        )
        writer = stk.XyzWriter()
        temp_string = writer.to_string(molecule)
        string += '\n'.join(temp_string.split('\n')[2:])
        string += '&END COORD\n\n'
        return string

    def _basis_set_section(self):
        return (
            '! Atom Basis Set\n'
            '! --------------\n'
            '&KIND C\n'
            '  BASIS_SET DZVP-MOLOPT-GTH\n'
            '  POTENTIAL GTH-PBE-q4\n'
            '&END KIND\n'
            '&KIND H\n'
            '  BASIS_SET DZVP-MOLOPT-GTH\n'
            '  POTENTIAL GTH-PBE-q1\n'
            '&END KIND\n'
            '&KIND N\n'
            '  BASIS_SET DZVP-MOLOPT-GTH\n'
            '  POTENTIAL GTH-PBE-q5\n'
            '&END KIND\n'
            '&KIND O\n'
            '  BASIS_SET DZVP-MOLOPT-GTH\n'
            '  POTENTIAL GTH-PBE-q6\n'
            '&END KIND\n'
            '&KIND Zn\n'
            '  BASIS_SET  DZVP-MOLOPT-SR-GTH\n'
            '  POTENTIAL  GTH-PBE-q12\n'
            '&END KIND\n'
            '&END SUBSYS\n'
            '&END FORCE_EVAL\n\n'
        )

    def _geomopt_section(self):
        raise NotImplementedError()

    def write_calculation_input(self):
        raise NotImplementedError()


class CP2KEnergy(CP2KCalculator):

    def _global_section(self):
        string = (
            '&GLOBAL\n'
            f'    PROJECT_NAME {self._job_name}\n'
            '    RUN_TYPE Energy\n'
            '    EXTENDED_FFT_LENGTHS\n'
            '    PRINT_LEVEL MEDIUM\n'
            '&END GLOBAL\n\n'
        )
        return string

    def write_calculation_input(
        self,
        output_directory,
        molecule,
        charge,
        guess,
        cutoff,
        rel_cutoff,
        solvent=None,
    ):
        string = self._global_section()
        string += self._force_section(
            charge=charge,
            cutoff=cutoff,
            rel_cutoff=rel_cutoff,
            guess=guess,
            solvent=solvent,
        )
        string += self._cell_section(molecule)
        string += self._coordinates_section(molecule)
        string += self._basis_set_section()

        input_file = os.path.join(
            output_directory, f'{self._job_name}.in'
        )
        with open(input_file, 'w') as f:
            f.write(string)


class CP2KOptimizer(CP2KCalculator):

    def _global_section(self):
        string = (
            '&GLOBAL\n'
            f'  PROJECT_NAME {self._job_name}\n'
            '  RUN_TYPE GEO_OPT\n'
            '  EXTENDED_FFT_LENGTHS\n'
            '&END GLOBAL\n\n'
        )
        return string

    def _geomopt_section(self):
        string = (
            '! Geometry Optimization\n'
            '! ---------------------\n'
            '&MOTION\n'
            '  &GEO_OPT\n'
            '    OPTIMIZER BFGS\n'
            '    MAX_ITER 2000\n'
            '    MAX_DR  [bohr] 0.003\n'
            '    &BFGS\n'
            '    &END\n'
            '  &END GEO_OPT\n'
            '&END MOTION\n\n'
        )
        return string

    def write_calculation_input(
        self,
        output_directory,
        molecule,
        charge,
        guess,
        cutoff,
        rel_cutoff,
        solvent=None,
    ):
        string = self._global_section()
        string += self._force_section(
            charge=charge,
            cutoff=cutoff,
            rel_cutoff=rel_cutoff,
            guess=guess,
            solvent=solvent,
        )
        string += self._cell_section(molecule)
        string += self._coordinates_section(molecule)
        string += self._basis_set_section()
        string += self._geomopt_section()

        input_file = os.path.join(
            output_directory, f'{self._job_name}.in'
        )
        with open(input_file, 'w') as f:
            f.write(string)


def write_sub_file(input_files, dft_directory):
    orca_bin_dir = '/apps/orca/4.2.1/bin/orca'

    runfile = os.path.join(f'{dft_directory}', 'run_orca.sh')

    runlines = ''.join([
        f"{orca_bin_dir} {infile} > {infile.replace('.in', '.out')}\n"
        for infile in input_files
    ])

    string = (
        '#PBS -N lig_spe\n'
        '#PBS -l walltime=72:00:00\n'
        '#PBS -l select=1:ncpus=32:mem=124gb\n\n'
        'module load orca/4.2.1\n\n'
        'cd $PBS_O_WORKDIR\n\n'
        f'{runlines}'
    )

    with open(runfile, 'w') as f:
        f.write(string)


def write_molecule_section(directory, base_name, struct, charge):

    multiplicity = 1
    xyzfile = f'{base_name}_init.xyz'

    string = f'* xyzfile {charge} {multiplicity} {xyzfile}\n'

    struct.write(os.path.join(directory, xyzfile))

    return string


def write_input_file(
    input_file,
    dft_directory,
    method,
    mol_file,
    charge,
):
    """
    Write ORCA single point energy of mol file.

    """

    if method == 'spe-pbe':
        top_line = (
            '! DFT SP RKS PBE0 def2-SVP D4 '
            'Grid6 NOFINALGRID SlowConv TightSCF\n\n'
        )
    elif method == 'opt-b97':
        top_line = (
            '! DFT COPT B97-3c '
            'Grid6 NOFINALGRID SlowConv TightSCF\n\n'
        )
    else:
        raise ValueError(
            'Method options:\n'
            'spe-pbe: ! DFT SP RKS PBE0 def2-SVP D4 Grid6 NOFINALGRID '
            'SlowConv TightSCF\n'
            'opt-b97: ! DFT COPT B97-3c Grid6 NOFINALGRID SlowConv '
            'TightSCF\n'
        )

    base_name = '_'+input_file.replace('.in', '')
    base_line = f'%base "{base_name}" \n\n'

    scf_section = (
        '%scf\n   MaxIter 2000\nend\n\n'
    )

    procs_section = (
        '%pal\n   nprocs 32\nend\n\n'
    )

    struct = stk.BuildingBlock.init_from_file(mol_file)
    mol_section = write_molecule_section(
        dft_directory, base_name, struct, charge
    )

    string = top_line
    string += base_line
    string += scf_section
    string += procs_section
    string += mol_section

    with open(f'{dft_directory}/{input_file}', 'w') as f:
        f.write(string)
