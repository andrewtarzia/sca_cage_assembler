"""Modules defining and building the Cage class and subclasses."""

import stk
from spindle_graph import M8L6KnotPrism


def m8l6_graph(building_blocks, vertex_alignments):
    topology_graph = stk.cage.M8L6Cube(
        building_blocks=building_blocks,
        vertex_alignments=vertex_alignments,
        num_processes=2,
    )

    return topology_graph


def m8l6knot_graph(building_blocks, vertex_alignments):
    topology_graph = M8L6KnotPrism(
        building_blocks=building_blocks,
        vertex_alignments=vertex_alignments,
        num_processes=2,
    )

    return topology_graph


def available_topologies(string):
    """Get stk function of desired topology."""
    topologies = {
        "m8l6face": m8l6_graph,
        "m8l6knot": m8l6knot_graph,
    }

    try:
        return topologies[string]
    except KeyError as e:
        msg = f"{string} not in {topologies.keys()}"
        raise KeyError(msg) from e


def metal_FFs(CN):
    """Define metal FF names for UFF4MOF.

    Key = Atomic number
    Value = UFF4MOF type
    CN = coordination number of metal.

    """
    # Default settings.
    dicts = {
        26: "Fe4+2",
        27: "Co4+2",
        28: "Ni4+2",  # No alternative available.
        30: "Zn4+2",  # No alternative available for 90 degrees.
        42: "Mo4f2",
        45: "Rh6+3",  # No alternative available.
        46: "Pd4+2",
        48: "Cd4f2",  # No alternative available.
        78: "Pt4+2",
    }

    if CN == 4:  # noqa: PLR2004
        pass
    elif CN == 6:  # noqa: PLR2004
        dicts[26] = "Fe6+2"
        dicts[27] = "Co6+2"

    return dicts
