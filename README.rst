:author: Andrew Tarzia

Overview
========

A Python module for generating and analysing metal-organic cages formed
from subcomponent self-assembly.

Code is for the paper at DOI: XX

If you have any issues contact me at ``andrew dot tarzia at gmail dot com``

Installation
============

To install this code, you can clone this repo. It is written to be used following the workflow below.

Dependancies:

The file ``environment.yml`` defines a Python 3.9 conda env used for this project.

The file ``env_set.py`` sets the directories used in this work based on my workstation (Ubuntu 18, 20 or 22). This will require editting for any other machine.

The external software used:

``CREST`` version ``2.9`` from https://github.com/crest-lab/crest

``Gulp`` version ``5.1`` from http://gulp.curtin.edu.au/gulp/

``Shape`` version ``2.1`` from http://www.ee.ub.edu/index.php?option=com_content&view=article&id=575:shape-available&catid=80:news&Itemid=466


Workflow used
=============

**Note** The file ``env_set.py`` sets the directories used in this work based on my workstation (Ubuntu 18, 20 or 22). This will require editting for any other machine.

``CODEDIR`` must be set by the user - where the code is.

``PROJDIR`` must be set by the user - where the user wants the data to be generated.

I attempt to highlight steps based on their part in the process (e.g. construction or analysis). Not all steps are necessary for what a user may want to do. Although, this project was so involved that many parts of the code are inter-related in an inconveniant to new usage. If so, please contact me for suggestions.

For visualisation.
------------------

plot shape measure:
    Directory to run in: ``PROJDIR``/shape_tests/

    Command:
    .. code-block::

        python CODEDIR/cage_structures/scripts/plot_shape_measure_examples.py


For ligand construction and analysis.
------------------

complex building:
    Directory to run in: ``PROJDIR``/complex_library/

    Produces: complex building block structures.

    Command: ``python CODEDIR/cage_structures/scripts/build_complex_library.py CODEDIR/cage_structures/data/cube_complex_library.json ../ligand_library/``


ligand building:
    Directory to run in: ``PROJDIR``/ligand_library/

    Produces: ligand building block structures.

    Command: ``python CODEDIR/cage_structures/scripts/build_ligand_library.py CODEDIR/cage_structures/data/cube_ligand_library.json``

Ligand AR analysis:
    Directory to run in: ``PROJDIR``/ligand_library/

    Performs: analysis of the aspect ratio of the ligands before and after reaction.

    Command: ``python CODEDIR/cage_structures/scripts/calculate_ligand_ARs.py``


flex analysis:
    Directory to run in: ``PROJDIR``/flex_analysis/

    Performs: analysis of the conformer flex of all ligands

    Command: ``python CODEDIR/cage_structures/scripts/flexibility_analysis.py ../ CODEDIR/cage_structures/data/cube_ligand_library.json``


For face construction.
------------------

face analysis:
    Directory to run in: ``PROJDIR``/face_analysis/

    Produces: face models and analyses them.

    Command: ``python CODEDIR/cage_structures/scripts/face_analysis.py ../ cl1 manual_complex``



For cage construction.
------------------

cage building:
    Directory to run in: ``PROJDIR``/cage_library/

    Produces: cage structures, optimised.

    Performs: analysis and saves them to a csv.

    Command: ``python CODEDIR/cage_structures/scripts/build_cube_library.py CODEDIR/cage_structures/data/cube_ligand_library.json CODEDIR/cage_structures/data/cube_complex_library.json CODEDIR/cage_structures/data/cube_library.json ../ligand_library/ ../complex_library/ f CODEDIR/cage_structures/data/cube_expt_library.json``

report on constructions:
    Directory to run in: ``PROJDIR``/cage_library/

    Produces: pdbs of optimised structures and a text file with report

    Command: ``python CODEDIR/cage_structures/scripts/report_on_construction.py``


For analysing and comparing to crystal structures (if available)
------------------

crystal structure analysis:
    Directory to run in: ``PROJDIR``/xray_structures/analysis/

    Performs: analysis of crystal structures using same methods as computational models.

    Command: ``python CODEDIR/cage_structures/scripts/analyse_crystal_structures.py CODEDIR/cage_structures/data/cube_complex_library.json CODEDIR/cage_structures/data/cube_library.json ../../ligand_library/ ../../cage_library/  CODEDIR/cage_structures/data/cube_expt_library.json``

align xray and generated structures:
    Directory to run in: ``PROJDIR``/alignment/

    Produces: many possible pairs of structures that are aligned, covering the multiple input rotations.

    Command: ``python CODEDIR/cage_structures/scripts/align_structures.py  CODEDIR/cage_structures/data/cube_complex_library.json CODEDIR/cage_structures/data/cube_library.json ../../../cage_library/  CODEDIR/cage_structures/data/cube_expt_library.json``

map pores of aligned xray and generated structures:
    Directory to run in: ``PROJDIR``/alignment/

    Produces: _pore.xyz and _host.xyz for each ccrystal structure

    Command: ``python CODEDIR/cage_structures/scripts/poremapping.py CODEDIR/cage_structures/data/cube_expt_library.json``


For setting up and performining DFT.
------------------

setup convergence tests:
    Directory to run in: ``PROJDIR``/cage_library/

    Produces: directory (set_dft_run) with input files for DFT energy evaluation as a function of parameters.

    Command: ``python CODEDIR/cage_structures/scripts/setup_convergence_tests.py conv_tests_dft ./ f``

evaluate convergence tests:
    Directory to run in: ``PROJDIR``/cage_library/

    Produces: plots of rel. energy in kJmol-1 vs cutoff or rel_cutoff

    Command: ``python CODEDIR/cage_structures/scripts/evaluate_convergence_tests.py conv_tests_dft``

setup set opt:
    Directory to run in: ``PROJDIR``/cage_library/

    Produces: directory (set_dft_run) with input files for CP2K DFT run.

    Command: ``python CODEDIR/cage_structures/scripts/setup_set_opt.py set_dft_run ./ cl1_quad2_12 f``

extract set opt:
    Directory to run in: ``PROJDIR``/cage_library/

    Produces: cage structures with _optdft.mol suffix

    Command: ``python CODEDIR/cage_structures/scripts/extract_set_opt.py ./set_dft_run ./ cl1_quad2_12``


For analysis/plotting.
------------------

plot flex measures:
    Directory to run in: ``PROJDIR``/flex_analysis/

    Produces: flex_dists.pdf and flex_comp.pdf and flex_energy.pdf

    Command: ``python CODEDIR/cage_structures/scripts/plot_flex_measures.py``

plot face measure examples:
    Directory to run in: ``PROJDIR``/cage_library/

    Produces: plots of simple models relationships between AR and face stability

    Command: ``python CODEDIR/cage_structures/scripts/plot_face_measure_examples.py``


plot categorisation:
    Directory to run in: ``PROJDIR``/cage_library/

    Produces:: categorical_*.pdf

    Command: ``python CODEDIR/cage_structures/scripts/plot_categorisation.py ../xray_structures/analysis/all_xray_csv_data.csv``


plot parities:
    Directory to run in: ``PROJDIR``/cage_library/

    Produces:: parities_*.pdf

    Command: ``python CODEDIR/cage_structures/scripts/plot_parities.py ../xray_structures/analysis/all_xray_csv_data.csv CODEDIR/cage_structures/data/cube_expt_library.json``

plot cube vs properties:
    Directory to run in: ``PROJDIR``/cage_library/

    Produces: shape_vs_energies.pdf and shape_vs_int_angle.pdf

    Performs: comparison of shape measure (cube likeness) with formation and strain energy

    Command: ``python CODEDIR/cage_structures/scripts/plot_cube_vs_properties.py``


plot lse vs fe:
    Directory to run in: ``PROJDIR``/cage_library/

    Produces: lse_sum_vs_fe.pdf and lse_sum_vs_fe_z.pdf

    Command: ``python CODEDIR/cage_structures/scripts/plot_lse_vs_fe.py``


plot set distributions:
    Directory to run in: ``PROJDIR``/cage_library/

    Produces: distribution_*pdf and set_energies_xtb/dft.pdf plots

    Command: ``python CODEDIR/cage_structures/scripts/plot_set_distributions.py``

plot symm distributions:
    Directory to run in: ``PROJDIR``/cage_library/

    Produces: sym_distribution_*.pdf figures

    Command: ``python CODEDIR/cage_structures/scripts/plot_symm_distributions.py``

decision tree:
    Directory to run in: ``PROJDIR``/cage_library/

    Produces: a decision tree plot — decision_tree.pdf

    Command: ``python CODEDIR/cage_structures/scripts/decision_tree.py``


plot znzn distributions:
    Directory to run in: ``PROJDIR``/cage_library/

    Produces: plots of zn-Zn distances for constructed and crystal structures.

    Command: ``python CODEDIR/cage_structures/scripts/plot_znzn_distributions.py ../xray_structures/analysis CODEDIR/cage_structures/data/cube_expt_library.json``

plot ligand properties:
    Directory to run in: ``PROJDIR``/cage_library/

    Produces: all_ligand_MM_vs_AR.pdf and all_ligand_properties.pdf

    Command: ``python CODEDIR/cage_structures/scripts/plot_ligand_properties.py CODEDIR/cage_structures/data/cube_expt_library.json``

plot td tl parity:
    Directory to run in: ``PROJDIR``/cage_library/

    Produces: td_tl parity plots.

    Command: ``python CODEDIR/cage_structures/scripts/plot_td_tl_parity.py``


Acknowledgements
================

I developed this code when I was working in the Jelfs group,
http://www.jelfs-group.org/.