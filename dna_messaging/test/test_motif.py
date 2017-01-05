#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dna_messaging.motif import *


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def test_motif_enumerations():
    assert motif_enumerations(['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT'], 3, 1) == 'ATA ATT GTT TTT'
