#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import copy
import sys

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


lines = sys.stdin.read().splitlines()
k = int(lines[0])
print(median_string(lines[1:], k))
