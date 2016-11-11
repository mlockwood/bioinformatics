__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


import sys


def reverse_complement(string):
    lookup = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    text = ''
    i = -1
    while i >= -len(string):
        text += lookup[string[i]]
        i -= 1
    return text


# lines = sys.stdin.read().splitlines()
# print(reverse_complement(*lines))
