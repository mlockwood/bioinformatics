#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dna_messaging.reverse_complement import *

__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def test_reverse_complement():
    assert reverse_complement('AAAACCCGGT') == 'ACCGGGTTTT'
