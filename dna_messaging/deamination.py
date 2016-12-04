#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def find_skew(genome):
    """
    Find the skew of a genome.
    :param genome: genome string
    :return: space delimited skew string
    """
    skew = [0]
    for base in genome:
        diff = 1 if base == 'G' else (-1 if base == 'C' else 0)
        skew.append(skew[-1] + diff)
    return ' '.join(str(i) for i in skew)


lines = sys.stdin.read().splitlines()
print(find_skew(lines))