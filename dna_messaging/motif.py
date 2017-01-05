#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import copy
import sys

from dna_messaging.constants import PROFILE_INDEX
from dna_messaging.generic import get_all_mismatched_kmers

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


def most_probable_in_profile(genome, k, profile_matrix):
    """
    Find the most probably motif given a profile matrix and a length of
    k within the genome.
    :param genome: genome string
    :param k: length of kmers/motif
    :param profile_matrix: profile probabilities for each base
    :return: most probably motif
    """
    max_score = (0, None)
    profile = get_profile_dict(profile_matrix)

    # Explore the genome by kmer of length k for each index until the end
    i = 0
    while i < len(genome) - int(k) + 1:

        # Consider the kmer by base, stop if the current sum is less than the max score
        j = 0
        score = 1
        while j < int(k):

            # Multiple the score by the profile probability of the current j index and observed base
            score *= profile[j][genome[i+j]]

            # If score probability is less than the max_score, stop processing kmer
            if score < max_score[0]:
                j = k
            j += 1

        # If score exceeded the max_score replace with the current kmer and its probability
        if score > max_score[0]:
            max_score = (score, genome[i:i+int(k)])
        i += 1

    return max_score[1]


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


lines = sys.stdin.read().splitlines()
genome = lines[0]
k = int(lines[1])
print(most_probable_in_profile(genome, k, lines[2:]))
