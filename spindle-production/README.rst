:author: Andrew Tarzia

Overview
========

Structure generation code for the paper at DOI: TBD

If you have any issues contact me at ``andrew dot tarzia at gmail dot com``

Workflow used
=============

:CODEDIR: must be set by the user - where the code is.
:PROJDIR: must be set by the user - where the user wants the data to be generated.

The spindles project reuses data from the 2022 paper, hence many steps can be
skipped.

Additionally, the code now automatically finds the data directory to provide
the library json files.


For extracting topologies.
--------------------------

extract topology:
    Produces: the topology graph, which can be written into the code as in
    ``spindle_graph.py``.

    .. code-block::

        python CODEDIR/spindle-production/extract_topology.py


For cage construction.
----------------------

cage building:
    Produces: cage structures, optimised.

    Performs: analysis and saves them to a csv.

    .. code-block::

        python CODEDIR/spindle-production/build_cube_library.py


For analysis/plotting.
----------------------

plot flex measures:
    Produces a series of figures comparing relative energies across previous
    work and this study.

    .. code-block::

        python CODEDIR/spindle-production/analyse_library.py
