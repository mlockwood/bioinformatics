__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


from comparing_genes.longest_common_substring import *


def test_solve_LCS():
    assert solve_LCS('AACCTTGG', 'ACACTGTGA') == 'AACTGG' or 'ACCTGG'
    assert solve_LCS('ABCXYZAYB', 'XYZABCB') == 'XYZAB'