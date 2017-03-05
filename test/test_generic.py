#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from generic import *

__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def test_reverse_complement():
    assert reverse_complement('AAAACCCGGT') == 'ACCGGGTTTT'


def test_get_all_kmers():
    assert set(get_all_kmers(2)) == {'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC',
                                     'TG', 'TT'}
    assert len(get_all_kmers(5)) == 1024


def test_hamming_distance():
    assert hamming_distance('GGGCCGTTGGT', 'GGACCGTTGAC') == 3


def test_get_all_mismatched_kmers():
    assert set(get_all_mismatched_kmers('AA', 1)) == {'AA', 'AC', 'AG', 'AT', 'CA', 'GA', 'TA'}


def test_lines_to_graph_dict():
    assert lines_to_graph_dict(['0 -> 1,2', '1 -> 0', '2 -> 0']) == {
        '0': {'1': True, '2': True},
        '1': {'0': True},
        '2': {'0': True}
    }
