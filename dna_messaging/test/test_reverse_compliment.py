__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


from dna_messaging.reverse_complement import *


def test_reverse_complement():
    assert reverse_complement('AAAACCCGGT') == 'ACCGGGTTTT'
