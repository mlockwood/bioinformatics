#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from genome_comparing.longest_common_substring import *


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


class TestSolveLCS:

    def test_same_length(self):
        assert solve_LCS('AACCTTGG', 'ACACTGTGA') == 'AACTGG' or 'ACCTGG'

    def test_different_length(self):
        assert solve_LCS('ABCXYZAYB', 'XYZABCB') == 'XYZAB'

    def test_global_alignment(self):
        assert solve_LCS('PLEASANTLY', 'MEANLY', method='global', sigma=5, mu=BLOSUM62
                         ) == (8, ('PLEASANTLY', '-MEA--N-LY')) or (8, ('PLEASANTLY', '-ME--AN-LY'))

    def test_local_alignment(self):
        assert solve_LCS('MEANLY', 'PENALTY', method='local', sigma=5, mu=PAM250
                         ) == (15, ('EANL-Y', 'ENALTY'))

    def test_fitting_alignment(self):
        assert solve_LCS('GTAGGCTTAAGGTTA', 'TAGATA', method='fitting', sigma=1, mu=1
            ) == (2, ('TAGGCTTA', 'TAGA--TA')) or (2, ('TAGGCTTA', 'TAG--ATA')) or (2, ('TAGGCTTA', 'TA-G-ATA'))

    def test_edit_distance(self):
        assert solve_LCS('PLEASANTLY', 'MEANLY', method='edit_distance') == 5


def test_DAG():
    graph = DAG.build_graph(['0->1:7', '0->2:4', '2->3:2', '1->4:1', '3->4:3'])
    assert(DAG(graph, '0', '4')).longest_path_in_DAG() == (9, ['0', '2', '3', '4'])
