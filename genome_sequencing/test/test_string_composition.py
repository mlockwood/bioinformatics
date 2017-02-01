#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from generic import *
from genome_sequencing.string_composition import *


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def test_string_composition():
    assert string_composition(5, 'CAATCCAAC') == 'AATCC\nATCCA\nCAATC\nCCAAC\nTCCAA'


def test_string_reconstruction_pre_ordered():
    assert string_reconstruction_pre_ordered('ACCGA\nCCGAA\nCGAAG\nGAAGC\nAAGCT') == 'ACCGAAGCT'


def test_get_prefix():
    assert get_prefix('TAA') == 'TA'


def test_get_suffix():
    assert get_suffix('TAA') == 'AA'


def test_get_overlap_graph():
    assert get_overlap_graph('ATGCG\nGCATG\nCATGC\nAGGCA\nGGCAT'
                             ) == 'AGGCA -> GGCAT\nCATGC -> ATGCG\nGCATG -> CATGC\nGGCAT -> GCATG\n'


def test_de_bruijn_graph_from_string():
    assert de_bruijn_graph_from_string(4, 'AAGATTCTCTAAGA') == (
        'AAG -> AGA,AGA\nAGA -> GAT\nATT -> TTC\nCTA -> TAA\nCTC -> TCT\nGAT -> ATT\nTAA -> AAG\nTCT -> CTA,CTC\n' +
        'TTC -> TCT\n'
    )


def test_de_bruijn_graph_by_composition():
    assert de_bruijn_graph_by_composition('GAGG\nCAGG\nGGGG\nGGGA\nCAGG\nAGGG\nGGAG\n'
                                          ) == 'AGG -> GGG\nCAG -> AGG,AGG\nGAG -> AGG\nGGA -> GAG\nGGG -> GGA,GGG\n'


def test_eulerian_cycle():
    """
    This can produce many different acceptable outputs so the test may fail but not actually prove an error in code.
    """
    assert eulerian_cycle(lines_to_graph_dict([
        '0 -> 3',
        '1 -> 0',
        '2 -> 1,6',
        '3 -> 2',
        '4 -> 2',
        '5 -> 4',
        '6 -> 5,8',
        '7 -> 9',
        '8 -> 7',
        '9 -> 6'
    ])) == ('6->8->7->9->6->5->4->2->1->0->3->2->6' or '9->6->5->4->2->1->0->3->2->6->8->7->9')
