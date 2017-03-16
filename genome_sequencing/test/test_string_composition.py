#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest

from genome_sequencing.string_composition import *


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


skip_random = pytest.mark.skipif(
    not pytest.config.getoption("--random"),
    reason="need --random option to run"
)


def test_string_composition():
    assert string_composition(5, 'CAATCCAAC') == 'AATCC\nATCCA\nCAATC\nCCAAC\nTCCAA'


def test_string_reconstruction_pre_ordered():
    assert string_reconstruction_pre_ordered(['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC',' AAGCT']) == 'ACCGAAGCT'


def test_get_prefix():
    assert get_prefix('TAA') == 'TA'


def test_get_suffix():
    assert get_suffix('TAA') == 'AA'


def test_get_overlap_graph():
    assert print_overlap_graph(get_overlap_graph(['ATGCG', 'GCATG', 'CATGC', 'AGGCA', 'GGCAT'])
                               ) == 'AGGCA -> GGCAT\nCATGC -> ATGCG\nGCATG -> CATGC\nGGCAT -> GCATG\n'


def test_de_bruijn_graph_from_string():
    assert print_de_bruijn_graph(de_bruijn_graph_from_string(4, 'AAGATTCTCTAAGA')) == (
        'AAG -> AGA,AGA\nAGA -> GAT\nATT -> TTC\nCTA -> TAA\nCTC -> TCT\nGAT -> ATT\nTAA -> AAG\nTCT -> CTA,CTC\n' +
        'TTC -> TCT\n'
    )


def test_de_bruijn_graph_by_composition():
    assert print_de_bruijn_graph(
        de_bruijn_graph_by_composition(['GAGG', 'CAGG', 'GGGG', 'GGGA', 'CAGG', 'AGGG', 'GGAG'])
    ) == 'AGG -> GGG\nCAG -> AGG,AGG\nGAG -> AGG\nGGA -> GAG\nGGG -> GGA,GGG\n'


class TestEulerianPath:

    @skip_random
    def test_eulerian_cycle(self):
        """
        This can produce many different acceptable outputs so the test may fail but not actually prove an error in code.
        """
        assert eulerian_path(lines_to_graph_dict([
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

    def test_eulerian_path(self):
        assert eulerian_path(lines_to_graph_dict([
            '0 -> 2',
            '1 -> 3',
            '2 -> 1',
            '3 -> 0,4',
            '6 -> 3,7',
            '7 -> 8',
            '8 -> 9',
            '9 -> 6'
        ]), nearly=True) == '6->7->8->9->6->3->0->2->1->3->4'


def test_resolve_overlaps():
    assert resolve_overlaps('AA->AT->TT') == 'AATT'


def test_genome_reconstruction():
    assert genome_reconstruction(['CTTA', 'ACCA', 'TACC', 'GGCT', 'GCTT', 'TTAC']) == 'GGCTTACCA'


def test_universal_circular_string():
    assert universal_circular_string(2) == '0011' or '0110' or '1001' or '1100'


def test_string_spelled_by_gapped_patterns():
    assert strings_spelled_by_gapped_patterns([
        'GACC|GCGC',
        'ACCG|CGCC',
        'CCGA|GCCG',
        'CGAG|CCGG',
        'GAGC|CGGA'
    ], 4, 2) == 'GACCGAGCGCCGGA'


def test_read_pairs_path():
    assert read_pairs_path([
        'GAGA|TTGA',
        'TCGT|GATG',
        'CGTG|ATGT',
        'TGGT|TGAG',
        'GTGA|TGTT',
        'GTGG|GTGA',
        'TGAG|GTTG',
        'GGTC|GAGA',
        'GTCG|AGAT'
    ], 4, 2) == 'GTGGTCGTGAGATGTTGA'


def test_maximal_non_branching_paths():
    assert maximal_non_branching_paths(lines_to_graph_dict([
        '1 -> 2',
        '2 -> 3',
        '3 -> 4, 5',
        '6 -> 7',
        '7 -> 6'
    ])) == [(1, 2, 3), (3, 4), (3, 5), (6, 7, 6)] or [(1, 2, 3), (3, 4), (3, 5), (7, 6, 7)]