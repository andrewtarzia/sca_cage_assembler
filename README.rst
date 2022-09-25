:author: Andrew Tarzia

Overview
========

A Python module for generating and analysing metal-organic cages formed
from subcomponent self-assembly.

Code is for the paper at DOI: XX

If you have any issues contact me at `andrew dot tarzia at gmail dot com`

Installation
============

To install this code, you can clone this repo. It is written to be used following the workflow below.

Dependancies:

The file `environment.yml` defines a Python 3.9 conda env used for this project.

The file `env_set.py` sets the directories used in this work based on my workstation (Ubuntu 18, 20 or 22). This will require editting for any other machine.

The external software used:

`CREST` version `2.9` from https://github.com/crest-lab/crest

`Gulp` version `5,1` from http://gulp.curtin.edu.au/gulp/

`Shape` version `2.1` from http://www.ee.ub.edu/index.php?option=com_content&view=article&id=575:shape-available&catid=80:news&Itemid=466


Workflow used
=============

**Note** The file `env_set.py` sets the directories used in this work based on my workstation (Ubuntu 18, 20 or 22). This will require editting for any other machine.

I attempt to highlight steps based on their part in the process (e.g. construction or analysis). Not all steps are necessary for what a user may want to do. Although, this project was so involved that many parts of the code are inter-related in an inconveniant to new usage. If so, please contact me for suggestions.


plot_shape_measure
	in shape_tests/
	python ~/thesource/cage_structures/scripts/plot_shape_measure_examples.py


ligand building:
	in ligand_library/
		produces ligand building block structures.
	python ~/thesource/cage_structures/scripts/build_ligand_library.py ~/thesource/cage_structures/data/cube_ligand_library.json


Ligand AR analysis:
	in ligand_library/
		analyses the aspect ratio of the ligands before and after reaction.
	python ~/thesource/cage_structures/scripts/calculate_ligand_ARs.py


flex analysis:
	in flex_analysis/
		analyses the conformer flex of all ligands
	python ~/thesource/cage_structures/scripts/flexibility_analysis.py ../ ~/thesource/cage_structures/data/cube_ligand_library.json


plot flex measures
	in flex_analysis/
		produces flex_dists.pdf and flex_comp.pdf and flex_energy.pdf
	python ~/thesource/cage_structures/scripts/plot_flex_measures.py

complex building
	in complex_library/
		produces complex building block structures.
	python ~/thesource/cage_structures/scripts/build_complex_library.py ~/thesource/cage_structures/data/cube_complex_library.json ../ligand_library/

face analysis:
	in face_analysis/
		produces face models and analyses them.
	python ~/thesource/cage_structures/scripts/face_analysis.py ../ cl1 manual_complex

plot_face_measure_examples.py
	in cage_library/
		produces plots of simple models relationships between AR and face stability
	 python ~/thesource/cage_structures/scripts/plot_face_measure_examples.py

cage building
	in cage_library/
		produces cage structures, optimised.
		performs analysis and saves them to a csv.
	python ~/thesource/cage_structures/scripts/build_cube_library.py ~/thesource/cage_structures/data/cube_ligand_library.json ~/thesource/cage_structures/data/cube_complex_library.json ~/thesource/cage_structures/data/cube_library.json ../ligand_library/ ../complex_library/ f ~/thesource/cage_structures/data/cube_expt_library.json

report on constructions
	in cage_library/
		produces pdbs of optimised structures and a text file with report
	python ~/thesource/cage_structures/scripts/report_on_construction.py

crystal structure analysis
	in xray_structures/analysis/
		analyses crystal structures using same methods as computational methods.
	python ~/thesource/cage_structures/scripts/analyse_crystal_structures.py ~/thesource/cage_structures/data/cube_complex_library.json ~/thesource/cage_structures/data/cube_library.json ../../ligand_library/ ../../cage_library/  ~/thesource/cage_structures/data/cube_expt_library.json

align xray and generated structures
	in alignment/
		produces many possible pairs of structures that are aligned, covering the multiple input rotations.
	python ~/thesource/cage_structures/scripts/align_structures.py  ~/thesource/cage_structures/data/cube_complex_library.json ~/thesource/cage_structures/data/cube_library.json ../../../cage_library/  ~/thesource/cage_structures/data/cube_expt_library.json


map pores of aligned xray and generated structures
	in alignment/
		produces _pore.xyz and _host.xyz for each ccrystal structure
	python ~/thesource/cage_structures/scripts/poremapping.py ~/thesource/cage_structures/data/cube_expt_library.json

plot_categorisation
	in cage_library/
		produces: categorical_*.pdf
	python ~/thesource/cage_structures/scripts/plot_categorisation.py ../xray_structures/analysis/all_xray_csv_data.csv


plot_parities
	in cage_library/
		produces: parities_*.pdf
	python ~/thesource/cage_structures/scripts/plot_parities.py ../xray_structures/analysis/all_xray_csv_data.csv ~/thesource/cage_structures/data/cube_expt_library.json

plot_cube_vs_properties
	in cage_library/
		produces shape_vs_energies.pdf and shape_vs_int_angle.pdf
		comparison of shape measure (cube likeness) with formation and strain energy
	python ~/thesource/cage_structures/scripts/plot_cube_vs_properties.py


plot_lse_vs_fe
	in cage_library/
		produces lse_sum_vs_fe.pdf and lse_sum_vs_fe_z.pdf
	python ~/thesource/cage_structures/scripts/plot_lse_vs_fe.py


plot_set_distributions
	in cage_library/
		produces distribution_*pdf and set_energies_xtb/dft.pdf plots
	python ~/thesource/cage_structures/scripts/plot_set_distributions.py

plot_symm_distributions
	in cage_library/
		sym_distribution_*.pdf figures
	python ~/thesource/cage_structures/scripts/plot_symm_distributions.py

decision tree
	in cage_library/
		produces a decision tree plot — decision_tree.pdf
	python ~/thesource/cage_structures/scripts/decision_tree.py


plot_znzn_distributions
	in cage_library/
		produces plots of zn-Zn distances for constructed and crystal structures.
	python ~/thesource/cage_structures/scripts/plot_znzn_distributions.py ../xray_structures/analysis ~/thesource/cage_structures/data/cube_expt_library.json

plot_ligand_properties
	in cage_library/
		produces all_ligand_MM_vs_AR.pdf and all_ligand_properties.pdf
	python ~/thesource/cage_structures/scripts/plot_ligand_properties.py ~/thesource/cage_structures/data/cube_expt_library.json

plot_td_tl_parity
	in cage_library/
		produces td_tl parity plots.
	python ~/thesource/cage_structures/scripts/plot_td_tl_parity.py

setup_convergence_tests
	in cage_library/
		produces directory (set_dft_run) with input files for DFT energy evaluation as a function of parameters.
	python ~/thesource/cage_structures/scripts/setup_convergence_tests.py conv_tests_dft ./ f

evaluate_convergence_tests
	in cage_library/
		produces plots of rel. energy in kJmol-1 vs cutoff or rel_cutoff
	python ~/thesource/cage_structures/scripts/evaluate_convergence_tests.py conv_tests_dft

setup_set_opt
	in cage_library/
		produces directory (set_dft_run) with input files for CP2K DFT run.
	python ~/thesource/cage_structures/scripts/setup_set_opt.py set_dft_run ./ cl1_quad2_12 f

extract_set_opt
	in cage_library/
		produces cage structures with _optdft.mol suffix
	python ~/thesource/cage_structures/scripts/extract_set_opt.py ./set_dft_run ./ cl1_quad2_12


Acknowledgements
================

I developed this code when I was working in the Jelfs group,
http://www.jelfs-group.org/.