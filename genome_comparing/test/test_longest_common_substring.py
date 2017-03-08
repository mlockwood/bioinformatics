__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


from genome_comparing.longest_common_substring import *


def test_solve_LCS():
    assert solve_LCS('AACCTTGG', 'ACACTGTGA') == 'AACTGG' or 'ACCTGG'
    assert solve_LCS('ABCXYZAYB', 'XYZABCB') == 'XYZAB'


def test_DAG():
    graph = DAG.build_graph(['0->1:7', '0->2:4', '2->3:2', '1->4:1', '3->4:3'])
    assert(DAG(graph, '0', '4')).longest_path_in_DAG() == (9, ['0', '2', '3', '4'])
