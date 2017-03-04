#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

from constants import *


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def translate_rna_to_amino_acids(genome):
    """
    Take a RNA string and convert it to its corresponding amino acid
    chain.
    :param genome: RNA string
    :return: amino acid string
    """
    i = 0
    aminos = ''
    while i < len(genome):
        aminos += GENETIC_CODE[genome[i:i+3]] if GENETIC_CODE[genome[i:i+3]] != '*' else ''
        i += 3
    return aminos


lines = sys.stdin.read().rstrip().splitlines()
print(translate_rna_to_amino_acids(lines[0]))
