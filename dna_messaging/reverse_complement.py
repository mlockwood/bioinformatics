#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def reverse_complement(string):
    lookup = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    text = ''
    i = -1
    while i >= -len(string):
        text += lookup[string[i]]
        i -= 1
    return text

