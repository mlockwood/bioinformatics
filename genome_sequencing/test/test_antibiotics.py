#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from genome_sequencing.antibiotics import *

__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def test_transcribe_dna_to_rna():
    assert transcribe_dna_to_rna('ATTCG') == 'AUUCG'


def test_transcribe_rna_to_dna():
    assert transcribe_rna_to_dna('AUUCG') == 'ATTCG'


def test_translate_rna_to_peptides():
    assert translate_rna_to_peptides('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA') == 'MAMAPRTEINSTRING'


def test_find_encoded_peptides():
    assert find_encoded_peptides('ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA', 'MA') == [
        'ATGGCC',
        'ATGGCC',
        'GGCCAT'
    ]


def test_count_cyclic_subpeptides():
    assert count_cyclic_subpeptides(31315) == 980597910


class TestFindPeptideSpectrum():

    def test_linear_spectrum(self):
        assert find_peptide_spectrum('NQEL') == [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484]

    def test_cyclic_spectrum(self):
        assert find_peptide_spectrum('NQEL', cyclic=True) == [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370,
                                                              371, 484]
