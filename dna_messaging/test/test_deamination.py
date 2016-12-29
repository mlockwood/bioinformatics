#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from dna_messaging.deamination import *

__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def test_find_skew():
    assert find_skew('CATGGGCATCGGCCATACGCC') == '0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2'


def test_find_min_skews():
    assert find_min_skews('TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT') == '11 24'
