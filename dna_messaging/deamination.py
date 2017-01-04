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


def find_min_skews(genome):
    """
    Find the minimum skew indices of a genome.
    :param genome: genome string
    :return: space delimited string of minimum skews
    """
    prev = 0
    cur_min = 0
    min_skews = []
    i = 0
    while i < len(genome):
        prev += 1 if genome[i] == 'G' else (-1 if genome[i] == 'C' else 0)
        if prev < cur_min:
            cur_min = prev
            min_skews = [i+1]
        elif prev == cur_min:
            min_skews.append(i+1)
        i += 1
    return ' '.join(str(i) for i in min_skews)


def find_max_skews(genome):
    """
    Find the maximum skew indices of a genome.
    :param genome: genome string
    :return: space delimited string of maximum skews
    """
    prev = 0
    cur_max = 0
    max_skews = []
    i = 0
    while i < len(genome):
        prev += 1 if genome[i] == 'G' else (-1 if genome[i] == 'C' else 0)
        if prev > cur_max:
            cur_max = prev
            max_skews = [i+1]
        elif prev == cur_max:
            max_skews.append(i+1)
        i += 1
    return ' '.join(str(i) for i in max_skews)


if __name__ == "__main__":
    # lines = sys.stdin.read().splitlines()
    # print(find_min_skews(*lines))
    print(find_max_skews('CATTCCAGTACTTCATGATGGCGTGAAGA'))