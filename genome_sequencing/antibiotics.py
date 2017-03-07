#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

from constants import *
from generic import *


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def transcribe_dna_to_rna(dna):
    """
    Take a DNA string and convert it to its corresponding RNA string.
    :param dna: DNA string
    :return: RNA string
    """
    rna = ''
    for base in dna:
        rna += 'U' if base == 'T' else base
    return rna


def transcribe_rna_to_dna(rna):
    """
    Take a RNA string and convert it to its corresponding DNA string.
    :param rna: RNA string
    :return: DNA string
    """
    dna = ''
    for base in rna:
        dna += 'T' if base == 'U' else base
    return dna


def translate_rna_to_peptides(rna):
    """
    Take a RNA string and convert it to its corresponding amino acid
    peptide chain.
    :param rna: RNA string
    :return: amino acid peptide string
    """
    i = 0
    peptides = ''
    while i < len(rna):
        peptides += GENETIC_CODE[rna[i:i + 3]] if GENETIC_CODE[rna[i:i + 3]] != '*' else ''
        i += 3
    return peptides


def find_encoded_peptides(dna, peptides):
    """
    Take a DNA string and a peptide and find all occurrences in DNA of
    the peptide
    :param dna: DNA string
    :param peptides: amino acid peptide string
    :return: all k-mers that encode the peptide string where k ==
    peptides * 3
    """
    results = []
    # Find the RNA transcription of the DNA string and its reverse complement
    rna_strings = (transcribe_dna_to_rna(dna), transcribe_dna_to_rna(reverse_complement(dna)))
    reverse = False
    for rna in rna_strings:
        i = 0
        while i <= len(rna) - (len(peptides) * 3):
            success = True
            # For each possible peptide sequence
            for x in range(len(peptides)):
                # Check that the next three bases translate to the appropriate peptide in the sequence
                if GENETIC_CODE[rna[i+x*3:i+x*3+3]] != peptides[x]:
                    success = False
            # If all translations match the peptide string; convert them from RNA to DNA
            if success:
                # If it was part of the reverse complement, find the reverse of the k-mer
                if reverse:
                    results.append(reverse_complement(transcribe_rna_to_dna(rna[i:i + (len(peptides) * 3)])))
                else:
                    results.append(transcribe_rna_to_dna(rna[i:i + (len(peptides) * 3)]))

            i += 1
        reverse = True

    return sorted(results)


def count_cyclic_subpeptides(peptide):
    """
    Find the count of cyclic subpeptides for a peptide.
    :param peptide: amino acid peptide chain
    :return: count of subpeptides
    """
    if isinstance(peptide, str):
        peptide = len(peptide)
    return peptide * (peptide - 1)


def spectrum_str_to_list(spectrum):
    """
    Take a spectrum string and convert it to a list of integers.
    :param spectrum: space delimited spectrum string
    :return: spectrum as a list of integers
    """
    return sorted([int(s) for s in spectrum.split()] if isinstance(spectrum, str) else spectrum)


def find_peptide_spectrum(peptide, cyclic=False):
    """
    Find the theoretical mass spectrum for a peptide. If cyclic is
    False the result will be a linear spectrum. If True it will be a
    cyclic spectrum.

    Note: the peptide may either be a string of alphabet single letter
    sequenced peptides or a list of integers representing masses.
    :param peptide: amino acid peptide chain
    :param cyclic: toggle between linear or cyclic
    :return: theoretical mass spectrum for the peptide
    """
    # Handle string instance
    if isinstance(peptide, str):
        peptide = [PEPTIDE_TO_MASS[p] for p in peptide]

    # Find the mass of the peptide at each position
    prefix_mass = {0: 0}
    i = 0
    while i < len(peptide):
        prefix_mass[i+1] = peptide[i] + prefix_mass[i]
        i += 1

    # Save the full peptide mass
    peptide_mass = prefix_mass[len(peptide)]

    # Build the spectrum from the prefix masses
    spectrum = [0]
    i = 0
    while i < len(peptide):
        j = i + 1
        while j < len(peptide) + 1:
            spectrum.append(prefix_mass[j] - prefix_mass[i])

            # If cyclic process the elements of the spectrum that span from end to beginning
            if cyclic:
                if i > 0 and j < len(peptide):
                    spectrum.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))
            j += 1
        i += 1

    return sorted(spectrum)


def cyclopeptide_sequencing(spectrum):
    """
    Find all peptides that fit the theoretical mass spectrum provided.
    :param spectrum: theoretical mass spectrum
    :return: all matching peptides
    """
    # Handle string input
    spectrum = spectrum_str_to_list(spectrum)

    # Set up data structures and process
    spectrum_lookup = set(spectrum)
    peptides = [[]]
    out = []
    while peptides:
        # Expand peptides
        expanded = []
        for peptide in peptides:
            for mass in MASS_TO_PEPTIDE:
                if mass in spectrum_lookup:
                    expanded.append(peptide + [mass])

        # Check peptides
        peptides = []
        for peptide in expanded:
            # Check if the peptide's mass is equal to the spectrum's
            if sum(peptide) == spectrum[-1]:
                # Verify that the entire peptide string matches
                if find_peptide_spectrum(''.join(MASS_TO_PEPTIDE[i][0] for i in peptide), cyclic=True) == spectrum:
                    out.append(peptide)

            # Otherwise test if the peptide's spectrum fits the provided spectrum, if it does allow it to be expanded
            else:
                match = True
                for mass in find_peptide_spectrum(''.join(MASS_TO_PEPTIDE[i][0] for i in peptide)):
                    if mass not in spectrum_lookup:
                        match = False
                if match:
                    peptides.append(peptide)

    return sorted(out)


def leaderboard_cyclopeptide_sequencing(spectrum, n, amino_acids=MASS_TO_PEPTIDE):
    """
    Find the peptide that most closely matches the spectrum using a
    branch-and-bound strategy that keeps N leaders at each iteration.
    This algorithm will keep peptides in that already satisfied the
    metric to bump others out of the leadership bracket.
    peptide mass == spectrum mass
    :param spectrum: experimental mass spectrum
    :param n: amount of leaders kept at each iteration
    :param amino_acids: list of allowed amino_acid masses
    :return: highest scoring cyclopeptide
    """
    # Handle string input
    spectrum = spectrum_str_to_list(spectrum)

    # Convert spectrum from list to dict
    spectrum_lookup = {}
    for mass in spectrum:
        spectrum_lookup[mass] = spectrum_lookup.get(mass, 0) + 1

    # Process iterations
    peptides = {tuple(): True}
    lead_peptide = {0: []}
    while peptides:
        # Expand peptides
        for peptide in list(peptides.keys()):
            for mass in amino_acids:
                # Ensure that no peptide can be added that would in sum exceed the spectrum's mass
                if sum(peptide) + mass <= spectrum[-1]:
                    peptides[(*peptide, mass)] = True
            del peptides[peptide]

        # Check peptides
        for peptide in list(peptides.keys()):
            # Check if the peptide's mass is equal to the spectrum's
            if sum(peptide) == spectrum[-1]:
                # Check if the score for the current peptide is greater than or equal to the lead score
                score = score_peptide(peptide, spectrum_lookup, cyclic=True)
                prev_score = next(iter(lead_peptide))
                if score > prev_score:
                    lead_peptide[score] = [peptide]
                    del lead_peptide[prev_score]
                if score == prev_score:
                    lead_peptide[prev_score] = lead_peptide.get(prev_score) + [peptide]
                peptides[peptide] = score

            # Score peptide
            else:
                peptides[peptide] = score_peptide(peptide, spectrum_lookup)

        # Select the N best peptides and delete the rest
        scores = sorted([(peptides[peptide], peptide) for peptide in peptides], reverse=True)
        n = select_scores_with_ties(scores, n)

        # Delete all of the other peptides
        for peptide_score in scores[n:]:
            del peptides[peptide_score[1]]

    return lead_peptide[next(iter(lead_peptide))]


def convolutional_cyclopeptide_sequencing(spectrum, m, n):
    """
    Use the M most frequent masses in the convolution of the spectrum
    to find the highest scoring peptide(s) given the leaderboard method
    using N.
    :param spectrum: experimental mass spectrum
    :param m: number of masses to keep
    :param n: peptides to keep during each leaderboard iteration
    :return: highest scoring peptide(s)
    """
    # Handle string instance
    spectrum = spectrum_str_to_list(spectrum)

    # Find the M most frequent masses in the convolutional spectrum
    convolution = find_spectral_convolution(spectrum)
    scores = sorted([(convolution[mass], mass) for mass in convolution if 57 <= mass <= 200], reverse=True)
    m = select_scores_with_ties(scores, m)
    aminos = dict((mass[1], True) for mass in scores[:m])

    # Run leaderboard
    return leaderboard_cyclopeptide_sequencing(spectrum, n, amino_acids=aminos)


def print_cyclopeptide_sequences(peptides):
    """
    Take lists of peptide masses and print them as space delimited
    entries of hyphen delimited masses such as mass-mass mass-mass.
    :param peptides: list of peptide mass lists
    :return: string of delimited entries
    """
    return ' '.join(['-'.join(str(i) for i in peptide) for peptide in peptides])


def cyclopeptide_scoring(peptide, spectrum):
    """
    Score a cyclic peptide against a spectrum.
    :param peptide: amino acid peptide chain
    :param spectrum: experimental mass spectrum
    :return: score
    """
    # Handle string input
    spectrum = spectrum_str_to_list(spectrum)

    # Convert spectrum from list to dict
    spectrum_lookup = {}
    for mass in spectrum:
        spectrum_lookup[mass] = spectrum_lookup.get(mass, 0) + 1

    return score_peptide(peptide, spectrum_lookup, cyclic=True)


def score_peptide(peptide, spectrum_lookup, cyclic=False):
    """
    Score a cyclic peptide against a spectrum's lookup dictionary.
    :param peptide: amino acid peptide chain
    :param spectrum_lookup: experimental mass spectrum as {mass: count}
    :param cyclic: toggle between linear and cyclic peptides
    :return: score
    """
    # Convert peptide to spectrum and then convert to dict
    peptide_lookup = {}
    for mass in find_peptide_spectrum(peptide, cyclic=cyclic):
        peptide_lookup[mass] = peptide_lookup.get(mass, 0) + 1

    # Score peptide
    score = 0
    for mass in peptide_lookup:
        if mass in spectrum_lookup:
            # If peptide count is >= spectrum add spectrum count
            if peptide_lookup[mass] >= spectrum_lookup[mass]:
                score += spectrum_lookup[mass]
            # Otherwise add peptide count
            else:
                score += peptide_lookup[mass]

    return score


def find_spectral_convolution(spectrum, zeros=False):
    """
    Find the convolution of masses for a spectrum.
    :param spectrum: experimental mass spectrum
    :param zeros: determination about whether to allow zeros
    :return: {convolution_mass: count}
    """
    # Handle string input
    spectrum = sorted(spectrum_str_to_list(spectrum), reverse=True)

    # Build convolution
    convolution = {}
    i = 0
    while i < len(spectrum):
        j = i + 1
        while j < len(spectrum):
            # Skip zeros if zeros is False
            if spectrum[i] != spectrum[j] or zeros:
                convolution[spectrum[i]-spectrum[j]] = convolution.get(spectrum[i]-spectrum[j], 0) + 1
            j += 1
        i += 1

    return convolution


def print_spectral_convolution(convolution):
    """
    Print a convolution accounting for multiplicity.
    :param convolution: a convolution dictionary
    :return: space delimited string of the convolution
    """
    out = []
    for key in sorted(convolution.keys()):
        [out.append(key) for x in range(convolution[key])]
    return ' '.join(str(s) for s in out)


#lines = sys.stdin.read().rstrip().splitlines()
