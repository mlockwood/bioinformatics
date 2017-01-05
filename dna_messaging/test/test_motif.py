#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dna_messaging.motif import *


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


PROFILE = ['0.2 0.2 0.3 0.2 0.3', '0.4 0.3 0.1 0.5 0.1', '0.3 0.3 0.5 0.2 0.4', '0.1 0.2 0.1 0.1 0.2']


def test_motif_enumerations():
    assert motif_enumerations(['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT'], 3, 1) == 'ATA ATT GTT TTT'


def test_median_string():
    assert median_string(['AAATTGACGCAT', 'GACGACCACGTT', 'CGTCAGCGCCTG', 'GCTGAGCACCGG', 'AGTTCGGGACAG'], 3) == 'GAC'


def test_most_probable_in_profile():
    assert most_probable_in_profile('ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT', 5, PROFILE) == 'CCGAG'


def test_get_profile_dict():
    assert get_profile_dict(PROFILE) == {
        0: {'A': 0.2, 'C': 0.4, 'G': 0.3, 'T': 0.1},
        1: {'A': 0.2, 'C': 0.3, 'G': 0.3, 'T': 0.2},
        2: {'A': 0.3, 'C': 0.1, 'G': 0.5, 'T': 0.1},
        3: {'A': 0.2, 'C': 0.5, 'G': 0.2, 'T': 0.1},
        4: {'A': 0.3, 'C': 0.1, 'G': 0.4, 'T': 0.2},
    }
