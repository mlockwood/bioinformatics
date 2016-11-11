__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


from class_I.pattern_count import *


def test_pattern_count():
    assert pattern_count('GCGCG', 'GCG') == 2


def test_frequent_kmers():
    assert frequent_kmers('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4) == 'CATG GCAT'


def test_pattern_match():
    assert pattern_match('ATAT', 'GATATATGCATATACTT') == '1, 3, 9'


def test_clump_finding():
    assert clump_finding('CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA', 5, 50, 4
                         ) == 'CGACA GAAGA'
