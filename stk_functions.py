#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions that are useful for stk usage

Author: Andrew Tarzia

Date Created: 18 Mar 2019
"""

import stk


def expected_window(topo):
    e_wind = {'dodec': 12, '4p6': 4, '4p62': 4,
          '8p12': 6,  '6p9': 5, '2p3': 3,
          '4p4': 6, '1p1': 3, '2p2': 4}
    return e_wind[topo]


def is_collapse(topo, avg_diff, max_window_diam, cavity_size, no_window):
    expected_wind = expected_window(topo)
    if expected_wind == no_window:
        alpha = 4 * avg_diff / (max_window_diam * expected_wind)
        if alpha < 0.035 and cavity_size > 1:
            # not collapsed
            return False
        else:
            # unknown
            return None
    else:
        # collapsed
        return True


def get_asymmetry(data):
    """Calculate assymetry as defined in GA paper (Berardo)

    The sum of all the windows' pair differences represents the asymmetry
    of the individual, Asymmetry parameter in eqn (1)
    """
    window_sizes = data['windows']['diameters']
    total = 0
    for i, a in enumerate(window_sizes):
        for j, b in enumerate(window_sizes[i:]):
            if i != j+i:
                diff = abs(a - b)
                total += diff
    return total
