#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def difference(s, t):
    """
    Take two binary vectors of a population and find the resulting
    difference value.
    :param s: binary vector
    :param t: binary vector
    :return: difference value
    """
    matches = set()
    i = 0
    while i < len(t):
        j = 0
        while j < len(t):
            if i != j and t[i] != t[j]:
                matches.add((i, j))
            j += 1
        i += 1
    denominator = len(matches)

    for match in list(matches):
        if s[match[0]] == s[match[1]]:
            matches.remove(match)

    return len(matches) / denominator


print(difference((1, 1, 0, 0, 0), (0, 0, 1, 1, 0)))
print(difference((0, 0, 1, 1, 0, 0, 1, 0), (1, 1, 0, 0, 1, 1, 1, 1)))
