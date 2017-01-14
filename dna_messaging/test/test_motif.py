#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dna_messaging.motif import *


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


GENOMES1 = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
GENOMES2 = [
    'AGGCGGCACATCATTATCGATAACGATTCGCCGCATTGCC',
    'ATCCGTCATCGAATAACTGACACCTGCTCTGGCACCGCTC',
    'AAGCGTCGGCGGTATAGCCAGATAGTGCCAATAATTTCCT',
    'AGTCGGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG',
    'AACCGGACGGCAACTACGGTTACAACGCAGCAAGAATATT',
    'AGGCGTCTGTTGTTGCTAACACCGTTAAGCGACGGCAACT',
    'AAGCGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG',
    'AATTGAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA'
]
GENOMES3 = [
    'GCACATCATTAAACGATTCGCCGCATTGCCTCGATAGGCG',
    'TCATAACTGACACCTGCTCTGGCACCGCTCATCCGTCGAA',
    'AAGCGGGTATAGCCAGATAGTGCCAATAATTTCCTTCGGC',
    'AGTCGGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG',
    'AACCGGACGGCAACTACGGTTACAACGCAGCAAGAATATT',
    'AGGCGTCTGTTGTTGCTAACACCGTTAAGCGACGGCAACT',
    'AAGCTTCCAACATCGTCTTGGCATCTCGGTGTGTGAGGCG',
    'AATTGAACATCTTACTCTTTTCGCTTTCAAAAAAAAGGCG'
]
GENOMES4 = [
    'GCACATCATTATCGATAACGATTCATTGCCAGGCGGCCGC',
    'TCATCGAATAACTGACACCTGCTCTGGCTCATCCGACCGC',
    'TCGGCGGTATAGCCAGATAGTGCCAATAATTTCCTAAGCG',
    'GTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTGAGTCG',
    'GACGGCAACTACGGTTACAACGCAGCAAGAATATTAACCG',
    'TCTGTTGTTGCTAACACCGTTAAGCGACGGCAACTAGGCG',
    'GCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTGAAGCG',
    'AAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAAAATTG'
]
MOTIF1 = {
    0: {0: 'A', 1: 'G', 2: 'G', 3: 'C'},
    1: {0: 'T', 1: 'G', 2: 'T', 3: 'T'},
    2: {0: 'G', 1: 'T', 2: 'G', 3: 'T'},
}
MOTIF2 = {
    0: {0: 'C', 1: 'C', 2: 'C', 3: 'C', 4: 'C'},
    1: {0: 'A', 1: 'A', 2: 'A', 3: 'A', 4: 'A'},
    2: {0: 'G', 1: 'G', 2: 'A', 3: 'A', 4: 'A'},
}
PROFILE = ['0.2 0.2 0.3 0.2 0.3', '0.4 0.3 0.1 0.5 0.1', '0.3 0.3 0.5 0.2 0.4', '0.1 0.2 0.1 0.1 0.2']


def test_motif_enumerations():
    assert motif_enumerations(['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT'], 3, 1) == 'ATA ATT GTT TTT'


def test_median_string():
    assert median_string(['AAATTGACGCAT', 'GACGACCACGTT', 'CGTCAGCGCCTG', 'GCTGAGCACCGG', 'AGTTCGGGACAG'], 3) == 'GAC'


def test_most_probable_in_profile():
    assert most_probable_in_profile('ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT', 5, get_profile_dict(PROFILE)
                                    ) == 'CCGAG'


def test_get_profile_dict():
    assert get_profile_dict(PROFILE) == {
        0: {'A': 0.2, 'C': 0.4, 'G': 0.3, 'T': 0.1},
        1: {'A': 0.2, 'C': 0.3, 'G': 0.3, 'T': 0.2},
        2: {'A': 0.3, 'C': 0.1, 'G': 0.5, 'T': 0.1},
        3: {'A': 0.2, 'C': 0.5, 'G': 0.2, 'T': 0.1},
        4: {'A': 0.3, 'C': 0.1, 'G': 0.4, 'T': 0.2},
    }


class TestGreedyMotifSearch:

    def test_without_laplace(self):
        assert greedy_motif_search(GENOMES1, 3) == ('CAG', 'CAG', 'CAA', 'CAA', 'CAA')

    def test_with_laplace(self):
        assert greedy_motif_search(GENOMES1, 3, True) == ('TTC', 'ATC', 'TTC', 'ATC', 'TTC')

    def test_off_by_one_error_at_start(self):
        assert greedy_motif_search(GENOMES2, 5, True) == ('AGGCG', 'ATCCG', 'AAGCG', 'AGTCG', 'AACCG', 'AGGCG', 'AGGCG',
                                                          'AGGCG')

    def test_off_by_one_error_at_end(self):
        assert greedy_motif_search(GENOMES3, 5, True) == ('AGGCG', 'TGGCA', 'AAGCG', 'AGGCA', 'CGGCA', 'AGGCG', 'AGGCG',
                                                          'AGGCG')

    def test_profile_tiebreaking(self):
        assert greedy_motif_search(GENOMES4, 5, True) == ('GGCGG', 'GGCTC', 'GGCGG', 'GGCAG', 'GACGG', 'GACGG', 'GGCGC',
                                                          'GGCGC')


def test_build_profile():
    assert build_profile(MOTIF1) == {
        0: {'A': 0.25, 'C': 0.25, 'G': 0.5},
        1: {'G': 0.25, 'T': 0.75},
        2: {'G': 0.5, 'T': 0.5},
    }


def test_score_matrix():
    assert score_matrix(MOTIF1) == 5


def test_convert_matrix_to_tuple():
    assert convert_matrix_to_tuple(MOTIF2) == ('CAG', 'CAG', 'CAA', 'CAA', 'CAA')
