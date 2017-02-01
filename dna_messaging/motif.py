#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import copy
import random
import sys

from numpy.random import choice

from dna_messaging.constants import BASES, PROFILE_INDEX
from generic import get_all_mismatched_kmers, hamming_distance

__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def motif_enumerations(genomes, k, d):
    """
    Find all kmers that occur within every genome in genomes with at
    most a hamming distance of d.
    :param genomes: list of genome strings
    :param k: length of k-mers
    :param d: maximum hamming distance
    :return: all kmers that occur in ever genome
    """
    kmers = {}  # records whether existing kmer matches are found in the current line
    patterns = {}  # lookup of observed patterns to kmers that can still satisfy being found in every line
    lookup = {}  # efficiency to recalculate strings at the end of each line by pointing which strings drop which kmers
    init = True  # toggle to differentiate initial action on first genome to actions on the following genomes

    for genome in genomes:

        # Explore the genome by pattern of length k for each index until the end
        i = 0
        while i < len(genome) - int(k) + 1:
            # If current pattern not in patterns find its matching kmers
            if genome[i:i+int(k)] not in patterns:
                # Find all potential mismatched kmers
                patterns[genome[i:i+int(k)]] = get_all_mismatched_kmers(genome[i:i+int(k)], int(d))

                # Setup lookup -> {kmer: {pattern: True}}
                if init:
                    # If this is the first line set all kmer options to the lookup
                    for kmer in patterns[genome[i:i+int(k)]]:
                        if kmer not in lookup:
                            lookup[kmer] = {}
                        lookup[kmer][genome[i:i+int(k)]] = True

                else:
                    # Only allow kmer options that already exist from previous lines
                    for kmer in list(patterns[genome[i:i+int(k)]].keys()):
                        if kmer in lookup:
                            lookup[kmer][genome[i:i+int(k)]] = True
                        # Remove kmers that have not appeared in every previous line from the pattern's lookup
                        else:
                            del patterns[genome[i:i+int(k)]][kmer]

            # For every kmer still matched to the pattern, set it to True for appearing in the line
            for kmer in patterns[genome[i:i+int(k)]]:
                kmers[kmer] = True

            i += 1

        # Handle end-of-line tasks
        for kmer in list(kmers.keys()):
            # If the kmer did not appear in this line remove it from all references
            if not kmers[kmer]:
                for pattern in lookup[kmer]:
                    del patterns[pattern][kmer]
                del lookup[kmer]
                del kmers[kmer]

            # Otherwise set it to False for processing on the next line
            else:
                kmers[kmer] = False

        init = False

    return ' '.join(sorted(kmers.keys()))


def median_string(genomes, k):
    """
    Find the k-mer with the minimum hamming distance score summed over
    all lines within the genomes.
    :param genomes: list of genome strings
    :param k: length of k-mers
    :return: motif k-mer with the lowest score the genome
    """
    kmers = {}  # records the score for all kmers
    patterns = {}  # lookup of observed patterns to kmers and their distance score
    lookup = {}  # store the minimum score for kmers in each line

    for genome in genomes:

        # Explore the genome by pattern of length k for each index until the end
        i = 0
        while i < len(genome) - int(k) + 1:
            # If current pattern not in patterns find its matching kmers
            if genome[i:i+int(k)] not in patterns:
                # Find all potential mismatched kmers
                patterns[genome[i:i+int(k)]] = get_all_mismatched_kmers(genome[i:i+int(k)], int(k))

            # If this is the first pattern for the genome line the lookup is equal to the pattern scores
            if not i:
                lookup = copy.deepcopy(patterns[genome[i:i+int(k)]])

            else:
                # Check to see if a new minimum applies to each kmer
                for kmer in patterns[genome[i:i+int(k)]]:
                    if patterns[genome[i:i+int(k)]][kmer] < lookup[kmer]:
                        lookup[kmer] = patterns[genome[i:i+int(k)]][kmer]

            i += 1

        # Handle end-of-line summing
        for kmer in lookup:
            kmers[kmer] = kmers.get(kmer, 0) + lookup[kmer]

    return min(kmers, key=kmers.get)


def get_profile_dict(profile_matrix):
    """
    Take a list of space delimited probabilities and convert it to a
    dictionary of {column: {base: probability}}.
    :param profile_matrix: list of space delimited probabilities
    :return: dictionary representation of profile matrix
    """
    profile = {}

    # Convert space delimited entries to actual lists
    matrix = []
    for row in profile_matrix:
        matrix.append([float(f) for f in row.split()])

    # Process by column
    j = 0
    while j < len(matrix[0]):

        # Add column to profile
        profile[j] = {}

        # Add each base and its probability
        i = 0
        while i < len(matrix):
            profile[j][PROFILE_INDEX[i]] = matrix[i][j]
            i += 1

        j += 1

    return profile


def greedy_motif_search(genomes, k, laplace=True):
    """
    Greedy algorithm solution to motif finding; find T motifs that
    produce the lowest score for the genome.
    :param genomes: set of genome strings
    :param k: length of kmers
    :param laplace: smoothing parameter for profiles
    :return: set of motifs with the lowest score
    """
    # Build initial best motifs from first kmers of each line
    best_motifs = dict((i, {}) for i in range(0, int(k)))
    i = 0
    while i < len(genomes):
        j = 0
        while j < int(k):
            best_motifs[j][i] = genomes[i][j]
            j += 1
        i += 1
    best_score = score_matrix(best_motifs)

    # Iterate through each kmer in the first genome string
    j = 0
    while j < len(genomes[0]) - int(k) + 1:
        # Set motif matrix to the first genome string's kmer
        motifs = dict((i, {0: genomes[0][j+i]}) for i in range(0, int(k)))

        # Iterate through each following genome string and update the motif and profile matrices
        i = 1
        while i < len(genomes):
            new_motif = most_probable_in_profile(genomes[i], int(k), build_profile(motifs, laplace))
            for x in range(0, int(k)):
                motifs[x][i] = new_motif[x]
            i += 1

        # If the score for the motifs is better than the best score set motifs as the best motifs
        if score_matrix(motifs) < best_score:
            best_motifs = motifs
            best_score = score_matrix(motifs)

        j += 1

    return convert_matrix_to_tuple(best_motifs)


def randomized_motif_search(genomes, k, laplace=True, trials=1000):
    """
    Randomized algorithm solution to motif finding; find T motifs that
    produce the lowest score for the genome.
    :param genomes: set of genome strings
    :param k: length of kmers
    :param laplace: smoothing parameter for profiles
    :param trials: number of times the random search will be completed
    :return: set of motifs with the lowest score
    """
    best_score = sys.maxsize
    best_motifs = {}

    # Process each trial
    for t in range(0, int(trials)):

        # Build initial best motifs from random kmers of each line
        motifs = dict((i, {}) for i in range(0, int(k)))
        i = 0
        while i < len(genomes):
            r = random.randint(0, len(genomes[i]) - int(k))
            j = 0
            while j < int(k):
                motifs[j][i] = genomes[i][r+j]
                j += 1
            i += 1

        # Iteratively find the most probable motifs given the previous profile until the score cannot be improved
        while motifs:
            profile = build_profile(motifs, laplace)

            # Iterate through each following genome string and update the motif
            i = 0
            new_motifs = dict((i, {}) for i in range(0, int(k)))
            while i < len(genomes):
                new_motif = most_probable_in_profile(genomes[i], int(k), profile)
                for x in range(0, int(k)):
                    new_motifs[x][i] = new_motif[x]
                i += 1

            # If the score for the motifs is better than the best score set motifs as the best motifs
            if score_matrix(new_motifs) < score_matrix(motifs):
                motifs = new_motifs

            # A better motif matrix was not found so terminate this trial
            else:
                # If this trial did better than all other previous trials record motifs and score
                if score_matrix(motifs) < best_score:
                    best_motifs = motifs
                    best_score = score_matrix(motifs)
                motifs = None

    return convert_matrix_to_tuple(best_motifs)


def randomized_gibbs_motif_search(genomes, k, gibbs=100, laplace=True, trials=1000):
    """
    Randomized algorithm solution to motif finding; find T motifs that
    produce the lowest score for the genome.
    :param genomes: set of genome strings
    :param k: length of kmers
    :param gibbs: number of times the gibbs search is run
    :param laplace: smoothing parameter for profiles
    :param trials: number of times the random search will be completed
    :return: set of motifs with the lowest score
    """
    best_score = sys.maxsize
    best_motifs = {}

    # Process each trial
    for t in range(0, int(trials)):

        # Build initial best motifs from random kmers of each line
        motifs = dict((i, {}) for i in range(0, int(k)))
        i = 0
        while i < len(genomes):
            r = random.randint(0, len(genomes[i]) - int(k))
            j = 0
            while j < int(k):
                motifs[j][i] = genomes[i][r+j]
                j += 1
            i += 1

        # Replace one motif at a time by random selection for gibbs number of times
        for g in range(0, int(gibbs)):
            # Select a random line for which to replace using an updated profile matrix
            r = random.randint(0, len(genomes) - 1)

            # Delete existing kmer motif of randomly selected line
            for col in motifs:
                del motifs[col][r]

            # Build a profile from the motifs and select a new best motif from the randomly selected line
            new_motif = random_probable_in_profile(genomes[r], int(k), build_profile(motifs, laplace))
            for x in range(0, int(k)):
                motifs[x][r] = new_motif[x]

            # If the score for the motifs is better than the best score set motifs as the best motifs
            if score_matrix(motifs) < best_score:
                best_motifs = copy.deepcopy(motifs)
                best_score = score_matrix(motifs)

    return convert_matrix_to_tuple(best_motifs)


def most_probable_in_profile(genome, k, profile):
    """
    Find the most probably motif given a profile matrix and a length of
    k within the genome.
    :param genome: genome string
    :param k: length of kmers/motif
    :param profile_matrix: profile probabilities for each base
    :return: most probable motif
    """
    max_prob = (0, genome[:int(k)])

    # Explore the genome by kmer of length k for each index until the end
    i = 0
    while i < len(genome) - int(k) + 1:

        # Consider the kmer by base, stop if the current sum is less than the max prob
        j = 0
        prob = 1
        while j < int(k):

            # Multiply the prob by the profile probability of the current j index and observed base
            try:
                prob *= profile[j][genome[i+j]]
            except KeyError:
                prob = 0

            # If prob probability is less than the max_prob, stop processing kmer
            if prob < max_prob[0]:
                j = k
            j += 1

        # If prob exceeded the max_prob replace with the current kmer and its probability
        if prob > max_prob[0]:
            max_prob = (prob, genome[i:i+int(k)])
        i += 1

    return max_prob[1]


def random_probable_in_profile(genome, k, profile):
    """
    Find a random motif given probabilities of a profile matrix and a
    length of k within the genome.
    :param genome: genome string
    :param k: length of kmers/motif
    :param profile_matrix: profile probabilities for each base
    :return: random probable motif
    """
    kmers = []
    probs = []
    total = 0.0

    # Explore the genome by kmer of length k for each index until the end
    i = 0
    while i < len(genome) - int(k) + 1:

        # Consider the kmer by base, stop if the current sum is less than the max prob
        j = 0
        prob = 1
        while j < int(k):

            # Multiply the prob by the profile probability of the current j index and observed base
            try:
                prob *= profile[j][genome[i+j]]
            except KeyError:
                prob = 0

            j += 1

        # Add kmer and probability to lists
        kmers.append(genome[i:i+k])
        probs.append(prob)
        total += prob
        i += 1

    # Normalize probabilities
    normalized_probs = []
    for prob in probs:
        normalized_probs.append(prob/total)

    return choice(kmers, p=normalized_probs)


def build_profile(motif_matrix, laplace=True):
    """
    Convert a motif_matrix to a profile matrix dict with probabilities.
    :param motif_matrix: a motif matrix in dictionary form
    :param laplace: smoothing parameter for profiles
    :return: profile matrix in dictionary form
    """
    profile = {}
    for col in motif_matrix:

        # Collect the distribution of bases for the column
        distr = dict((b, 1) for b in BASES) if laplace else {}
        for row in motif_matrix[col]:
            distr[motif_matrix[col][row]] = distr.get(motif_matrix[col][row], 0) + 1

        # Convert the distributions into probabilities
        for base in distr:
            distr[base] = distr.get(base, 0) / (len(motif_matrix[col]) + 4 if laplace else len(motif_matrix[col]))

        profile[col] = copy.deepcopy(distr)

    return profile


def score_matrix(motif_matrix):
    """
    Take a motif matrix in dictionary form {column: {row: base}} and
    output the score for the motif matrix.
    :param motif_matrix: a motif matrix in dictionary form
    :return: score for the motif matrix
    """
    score = 0
    for col in motif_matrix:

        # Collect the distribution of bases for the column
        distr = {}
        for row in motif_matrix[col]:
            distr[motif_matrix[col][row]] = distr.get(motif_matrix[col][row], 0) + 1

        # Score by taking the total rows less the count of the most frequent base
        score += (len(motif_matrix[col]) - distr[max(distr, key=distr.get)])
    return score


def convert_matrix_to_tuple(motif_matrix):
    """
    Take a motif matrix in dictionary form {column: {row: base}} and
    output a tuple of the motifs
    :param motif_matrix: a motif matrix in dictionary form
    :return: tuple of motifs
    """
    motifs = {}
    for col in sorted(motif_matrix.keys()):
        for row in motif_matrix[col]:
            motifs[row] = motifs.get(row, '') + motif_matrix[col][row]
    return tuple(v for k, v in sorted(motifs.items()))


def distance_between_pattern_and_strings(pattern, genomes):
    """
    Find the total distance between a pattern and each string in
    genomes. Take the minimum distance for each string. This is NOT
    used in median string because it is very inefficient.
    :param pattern: k-mer pattern
    :param genomes: set of genomes
    :return: distance
    """
    score = 0
    for genome in genomes:
        i = 0
        current = sys.maxsize
        while i <= len(genome) - len(pattern):
            if hamming_distance(pattern, genome[i:i+len(pattern)]) < current:
                current = hamming_distance(pattern, genome[i:i+len(pattern)])
            i += 1
        score += current
    return score
