#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from genome_comparing.manhattan_tourist import *


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def test_manhattan_tourist():
    reader = open('genome_comparing/test/manhattan_tourist_test_data.txt', 'r')
    lines = []
    for row in reader:
        lines.append(row.rstrip())
    assert process(lines) == 34
