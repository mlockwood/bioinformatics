__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


from class_I.reverse_complement import *


def test_reverse_complement():
    assert reverse_complement('AAAACCCGGT') == 'ACCGGGTTTT'
