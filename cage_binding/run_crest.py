#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run CREST (GFN) calculations.

Author: Andrew Tarzia

Date Created: 20 Mar 2019
"""
# run CREST
~/software/crest/crest pentacene_guest_1.xyz -xnam ~/software/xtb_190318/bin/xtb > pentacene_guest_1.output &
# with solvent
~/software/crest/crest pentacene_guest_1.xyz -g ch2cl2 -xnam ~/software/xtb_190318/bin/xtb > pentacene_guest_1.output &

# post CREST optimization
~/software/xtb_190318/bin/xtb XYZ -I xctrl --ohess > OUTPUT
# with solvent
~/software/xtb_190318/bin/xtb XYZ -I xctrl --ohess --gbsa > OUTPUT
