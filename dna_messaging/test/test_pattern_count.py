#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dna_messaging.pattern_count import *


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


STRING_1 = ('GAGCGCTGGGTTAACTCGCTACTTCCCGACGAGCGCTGTGGCGCAAATTGGCGATGAAACTGCAGAGAGAACTGGTCATCCAACTGAATTCTCCCCGCTATCGCAT'
'TTTGATGCGCGCCGCGTCGATT')
STRING_2 = ('CCAAATCCCCTCATGGCATGCATTCCCGCAGTATTTAATCCTTTCATTCTGCATATAAGTAGTGAAGGTATAGAAACCCGTTCAAGCCCGCAGCGGTAAAACCGAG'
'AACCATGATGAATGCACGGCGATTGCGCCATAATCCAAACA')
STRING_3 = ('CCGTCATCCGTCATCCTCGCCACGTTGGCATGCATTCCGTCATCCCGTCAGGCATACTTCTGCATATAAGTACAAACATCCGTCATGTCAAAGGGAGCCCGCAGCG'
'GTAAAACCGAGAACCATGATGAATGCACGGCGATTGC')


def test_pattern_count():
    assert pattern_count('GCGCG', 'GCG') == 2


def test_pattern_match():
    assert pattern_match('ATAT', 'GATATATGCATATACTT') == '1 3 9'


def test_clump_finding():
    assert clump_finding('CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA', 5, 50, 4
                         ) == 'CGACA GAAGA'


class TestApproximatePatternMatch:

    def test_sample(self):
        assert approximate_pattern_match(
            'ATTCTGGA',
            'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT',
            3
        ) == '6 7 26 27'

    def test_none_less_than_d(self):
        assert approximate_pattern_match('AAA', 'TTTTTTAAATTTTAAATTTTTT', 2) == '4 5 6 7 8 11 12 13 14 15'

    def test_only_less_than_d(self):
        assert approximate_pattern_match('TTT', 'AAAAAA', 3) == '0 1 2 3'

    def test_d_equal_none(self):
        assert approximate_pattern_match('CCA', 'CCACCT', 0) == '0'

    def test_first_index(self):
        assert approximate_pattern_match('GAGCGCTGG', STRING_1, 2) == '0 30 66'

    def test_last_indices(self):
        assert approximate_pattern_match('AATCCTTTCA', STRING_2, 3) == '3 36 74 137'

    def test_overlaps(self):
        assert approximate_pattern_match('CCGTCATCC', STRING_3, 3) == '0 7 36 44 48 72 79 112'


def test_count_approximate_pattern_match():
    assert count_approximate_pattern_matches('AAAAA', 'AACAAGCTGATAAACATTTAAAGAG', 1) == 4
    assert count_approximate_pattern_matches('AAAAA', 'AACAAGCTGATAAACATTTAAAGAG', 2) == 11


class TestGetKmerCounts:

    def test_zero_distance_counts(self):
        assert get_kmer_counts('AGGT', 2) == {'AA': 0, 'AC': 0, 'AG': 1, 'AT': 0, 'CA': 0, 'CC': 0, 'CG': 0, 'CT': 0,
                                              'GA': 0, 'GC': 0, 'GG': 1, 'GT': 1, 'TA': 0, 'TC': 0, 'TG': 0, 'TT': 0}

    def test_non_zero_distance_counts(self):
        assert get_kmer_counts('AGGT', 2, 1) == {'AA': 1, 'AC': 1, 'AG': 2, 'AT': 2, 'CA': 0, 'CC': 0, 'CG': 2, 'CT': 1,
                                                 'GA': 2, 'GC': 2, 'GG': 3, 'GT': 2, 'TA': 0, 'TC': 0, 'TG': 2, 'TT': 1}


def test_get_frequent_kmers():
    assert get_frequent_kmers('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4) == 'CATG GCAT'


def test_get_frequent_approximate_kmers():
    assert get_frequent_kmers('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 1) == 'ATGC ATGT GATG'


def test_get_frequent_approximate_and_reverse_kmers():
    assert get_frequent_kmers('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 1, rc=True) == 'ACAT ATGT'