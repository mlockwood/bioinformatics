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


class TestFindPeptideSpectrum:

    def test_linear_spectrum(self):
        assert find_peptide_spectrum('NQEL') == [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484]

    def test_cyclic_spectrum(self):
        assert find_peptide_spectrum('NQEL', cyclic=True) == [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370,
                                                              371, 484]


def test_cyclopeptide_sequencing():
    assert print_cyclopeptide_sequences(cyclopeptide_sequencing([0, 113, 128, 186, 241, 299, 314, 427])) == (
        '113-128-186 113-186-128 128-113-186 128-186-113 186-113-128 186-128-113'
    )


def test_leaderboard_cyclopeptide_sequencing():
    assert sorted(leaderboard_cyclopeptide_sequencing([0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460],
                                                      10)) == [
        (71, 129, 113, 147), (71, 147, 113, 129), (113, 129, 71, 147), (113, 147, 71, 129),
        (129, 71, 147, 113), (129, 113, 147, 71), (147, 71, 129, 113), (147, 113, 129, 71)
    ]


def test_convolutional_cyclopeptide_sequencing():
    assert (99, 71, 137, 57, 72, 57) in convolutional_cyclopeptide_sequencing(
        [57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285, 299, 307, 323, 356, 364, 394, 422, 493], 20, 60)


def test_cyclopeptide_scoring():
    assert cyclopeptide_scoring('NQEL', [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]) == 11


def test_find_spectral_convolution():
    assert print_spectral_convolution(find_spectral_convolution([0, 137, 186, 323])) == '49 137 137 186 186 323'
