#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Code to extract a topology from a molecule to make `knot_topology.py`.

Author: Andrew Tarzia

Date Created: 28 May 2021

"""

import stk
import stko
from rdkit.Chem import AllChem as rdkit


def extract_topo():
    struct = stk.BuildingBlock.init_from_file('initial_mol.mol')
    struct = struct.with_centroid([0, 0, 0])

    broken_bonds_by_id = []
    disconnectors = []
    smarts_to_search_for = '[#6r6X3]~[#7r5X3]~[#6r5X3H1]~[#6r5X3]'
    rdkit_mol = struct.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_mol)
    for atom_ids in rdkit_mol.GetSubstructMatches(
        query=rdkit.MolFromSmarts(smarts_to_search_for),
    ):
        bond_c = atom_ids[0]
        bond_n = atom_ids[1]
        broken_bonds_by_id.append(sorted((bond_c, bond_n)))
        disconnectors.extend((bond_c, bond_n))

    new_topology_graph = stko.TopologyExtractor()
    tg_info = new_topology_graph.extract_topology(
        molecule=struct,
        broken_bonds_by_id=broken_bonds_by_id,
        disconnectors=set(disconnectors),
    )
    struct.write('tg_cage.mol')
    tg_info.write('tg_info.pdb')
    return tg_info


def main():
    topo_info = extract_topo()

    v_pos = topo_info.get_vertex_positions()
    v_con = topo_info.get_connectivities()
    for i in v_pos:
        pos = v_pos[i]
        print(i, v_con[i], [round(i, 1) for i in pos])

    edges = topo_info.get_edge_pairs()
    for i, pair in enumerate(edges):
        print(
            'stk.Edge(\n'
            f'    id={i},\n'
            f'    vertex1=_vertex_prototypes[{pair[0]}],\n'
            f'    vertex2=_vertex_prototypes[{pair[1]}],\n'
            '),\n'
        )


if __name__ == '__main__':
    main()
