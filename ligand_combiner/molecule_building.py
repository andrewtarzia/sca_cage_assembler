#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build all molecules from DB of core, ligands and linkers.

Author: Andrew Tarzia

Date Created: 04 Apr 2019

"""

import sys
from os.path import join
from stk import rdkit_ETKDG
sys.path.insert(0, '/home/atarzia/thesource/')
from stk_functions import build_population, build_ABCBA, build_ABA


def main():
    if (not len(sys.argv) == 1):
        print("""
    Usage: molecule_builing.py
        """)
        sys.exit()
    else:
        # CIF = sys.argv[1]
        pass

    proj_dir = '/home/atarzia/projects/ligand_combiner/'
    core_dir = proj_dir + 'cores/'
    liga_dir = proj_dir + 'ligands/'
    link_dir = proj_dir + 'linkers/'
    mole_dir = proj_dir + 'molecules/'

    # build molecule populations
    core_pop = build_population(directory=core_dir, fgs=['bromine'])
    liga_pop = build_population(directory=liga_dir, fgs=['bromine'])
    link_pop = build_population(directory=link_dir, fgs=['bromine'])

    # build ABCBA polymer
    pop_ids = (0, 0, 0)  # core, ligand, linker
    large_part = build_ABCBA(core=core_pop[pop_ids[0]],
                             liga=liga_pop[pop_ids[1]],
                             link=link_pop[pop_ids[2]])
    # output as built
    prefix = core_pop[pop_ids[0]].name + '_'
    prefix += liga_pop[pop_ids[0]].name + '_'
    prefix += link_pop[pop_ids[0]].name
    json_file = prefix + '_ABCBA.json'
    large_part.dump(join(mole_dir, json_file))
    mol_file = prefix + '_ABCBA.mol'
    large_part.write(join(mole_dir, mol_file))
    # energy minimize polymer
    rdkit_ETKDG(large_part)
    # output energy minimized
    json_file = prefix + '_ABCBA_opt.json'
    large_part.dump(join(mole_dir, json_file))
    mol_file = prefix + '_ABCBA_opt.mol'
    large_part.write(join(mole_dir, mol_file))

    # build ABA polymer
    pop_ids = (0, 0, 0)  # core, ligand, linker
    small_part = build_ABA(core=core_pop[0],
                           liga=liga_pop[0])
    # output as built
    prefix = core_pop[pop_ids[0]].name + '_'
    prefix += liga_pop[pop_ids[0]].name
    json_file = prefix + '_ABA.json'
    small_part.dump(join(mole_dir, json_file))
    mol_file = prefix + '_ABA.mol'
    small_part.write(join(mole_dir, mol_file))
    # energy minimize polymer
    rdkit_ETKDG(small_part)
    # output energy minimized
    json_file = prefix + '_ABA_opt.json'
    small_part.dump(join(mole_dir, json_file))
    mol_file = prefix + '_ABA_opt.mol'
    small_part.write(join(mole_dir, mol_file))


if __name__ == "__main__":
    main()
