#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def manhattan_tourist(n, m, down, right):
    """
    Solve the manhattan tourist problem by finding the longest path to
    a grid where the only options are down and right.
    :param n: first dimension length
    :param m: second dimension length
    :param down: weights for moves down the matrix
    :param right: weights for moves right on the matrix
    :return: maximum weight for the entire matrix
    """
    path = {0: {0: 0}}

    i = 1
    while i <= n:
        path[i] = {0: (path[i-1][0] + down[i-1][0])}
        i += 1

    j = 1
    while j <= m:
        path[0][j] = (path[0][j-1] + right[0][j-1])
        j += 1

    i = 1
    while i <= n:
        j = 1
        while j <= m:
            path[i][j] = max([path[i-1][j] + down[i-1][j], path[i][j-1] + right[i][j-1]])
            j += 1
        i += 1

    return path[n][m]


def convert_matrix(matrix):
    """
    Convert a space and line delimited matrix to an indexed dictionary
    representation.
    :param matrix: space and line delimited matrix
    :return: indexed dictionary representation of the matrix
    """
    out = {}
    i = 0
    while i < len(matrix):
        out[i] = {}
        j = 0
        while j < len(matrix[i]):
            out[i][j] = int(matrix[i][j])
            j += 1
        i += 1
    return out


def process(lines):
    """
    Function to take the input variables and send them to the
    appropriate functions.
    :param lines: series of input variables
    :return: final result of the manhattan tourist problem given lines
    """
    n, m = re.split(' ', lines[0])
    n = int(n)
    m = int(m)

    down = []
    right = []
    down_toggle = True

    for line in lines[1:]:
        if re.search('-', line):
            down_toggle = False
        elif down_toggle:
            down.append(re.split(' ', line))
        else:
            right.append(re.split(' ', line))

    return manhattan_tourist(n, m, convert_matrix(down), convert_matrix(right))


# print(process(sys.stdin.read().splitlines()))
