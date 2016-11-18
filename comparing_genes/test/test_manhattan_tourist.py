__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


from comparing_genes.manhattan_tourist import *


def test_manhattan_tourist():
    reader = open('manhattan_tourist_test_data.txt', 'r')
    lines = []
    for row in reader:
        lines.append(row.rstrip())
    assert process(lines) == 34
