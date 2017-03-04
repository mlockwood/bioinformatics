#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from genome_sequencing.antibiotics import *

__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def test_translate_rna_to_amino_acids():
    assert translate_rna_to_amino_acids('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA') == 'MAMAPRTEINSTRING'
