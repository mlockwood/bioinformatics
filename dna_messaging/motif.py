#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
                patterns[genome[i:i+int(k)]] = get_all_mismatched_kmers(int(k), int(d), genome[i:i+int(k)])

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

        # Handle end-of-line functions
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

lines = sys.stdin.read().splitlines()
k, d = lines[0].split()
print(motif_enumerations(lines[1:], k, d))