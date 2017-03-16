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
        '0': {'1': 1, '2': 1},
        '1': {'0': 1},
        '2': {'0': 1}
    }
    assert lines_to_graph_dict(['0->1:7', '0->2:4', '2->3:2', '1->4:1', '3->4:3'], weighted=True) == {
        '0': {'1': 7, '2': 4},
        '1': {'4': 1},
        '2': {'3': 2},
        '3': {'4': 3}
    }


def test_select_scores_with_ties():
    assert select_scores_with_ties([(8, 'A'), (7, 'B'), (7, 'C'), (6, 'D')], 2) == 3


def test_list_of_dict_counts():
    assert list_of_dict_counts({1: 2, 3: 1, 2: 3}) == [1, 1, 2, 2, 2, 3]
