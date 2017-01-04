#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dna_messaging.generic import *

__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def test_get_all_kmers():
    assert set(get_all_kmers(2)) == {'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC',
                                     'TG', 'TT'}
    assert len(get_all_kmers(5)) == 1024


def test_hamming_distance():
    assert hamming_distance('GGGCCGTTGGT', 'GGACCGTTGAC') == 3
