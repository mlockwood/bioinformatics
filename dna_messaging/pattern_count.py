#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys

from dna_messaging.reverse_complement import reverse_complement


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


BASES = {'A', 'C', 'G', 'T'}


def get_all_kmers(k, kmer=''):
    """
    Generate a dictionary of all k-mers for k
    :param k: length of kmer
    :return: {kmer: 0}
    """
    kmers = {}
    if len(kmer) == k:
        return {kmer: 0}
    else:
        for base in BASES:
            kmers.update(get_all_kmers(k, '{}{}'.format(kmer, base)))
        return kmers


def pattern_count(genome, pattern):
    """
    Count the number of times a pattern appears in a genome genome.
    :param genome: genome string
    :param pattern: substring value to be searched for within the genome
    :return: count of pattern in genome string
    """
    count = 0
    i = 0
    while i <= len(genome) - len(pattern):
        if genome[i:i+len(pattern)] == pattern:
            count += 1
        i += 1
    return count


def pattern_match(substring, genome):
    """
    Return a series of indices in which a k-mer substring appears.
    :param substring: string of nucleotides that is being searched for
    :param genome: genome string
    :return: all indices of where the substring is found
    """
    index = []
    i = 0
    while i <= len(genome) - len(substring):
        if genome[i:i+len(substring)] == substring:
            index += [i]
        i += 1
    return ' '.join([str(i) for i in index])


def clump_finding(genome, k, L, t):
    """
    Identify all k-mer nucleotide strings that are found of at least
    frequency `t` within a window size (string distance) of `L`.
    :param genome: genome string
    :param k: length of nucleotide string
    :param L: window size
    :param t: minimum count threshold
    :return: all k-mers found within the window size at least `t` times
    """
    window = {}  # {k-mer: { index: True }}
    lookup = {}  # {index: k-mer}
    results = {}  # {k-mer: True}

    i = 0
    k = int(k)
    L = int(L)
    t = int(t)
    while i <= len(genome) - k:
        kmer = genome[i:i+k]

        # Add kmer to window
        if kmer not in window:
            window[kmer] = {}
        window[kmer][i] = True

        # Add index for kmer to lookup to manage window DS
        lookup[i] = kmer

        # Remove the kmer that became out of window
        if i >= L:
            del window[lookup[i-L]][i-L]

        # Add kmer to results if there are at least 't' occurrences within the range
        if len(window[kmer]) >= t:
            results[kmer] = True

        i += 1

    return ' '.join(sorted(list(results.keys())))


def hamming_distance(p, q):
    """
    Find the hamming distance between two strings.
    :param p: genome string
    :param q: genome string
    :return: hamming distance between the two input strings
    """
    i = 0
    distance = 0
    while i < len(p):
        if p[i] != q[i]:
            distance += 1
        i += 1
    return distance


def approximate_pattern_match(kmer, genome, d):
    """
    Find the indices where the kmer appears in the genome with a
    hamming distance up to `d`
    :param kmer: kmer string
    :param genome: genome string
    :param d: maximum distance value
    :return: space delimited list of indices
    """
    i = 0
    indices = []
    while i <= len(genome) - len(kmer):
        if hamming_distance(kmer, genome[i:i+len(kmer)]) <= int(d):
            indices.append(i)
        i += 1
    return ' '.join(str(i) for i in indices)


def count_approximate_pattern_matches(kmer, genome, d):
    """
    Find the count of how often the kmer appears in the genome with a
    hamming distance up to `d`
    :param kmer: kmer string
    :param genome: genome string
    :param d: maximum distance value
    :return: count
    """
    return len(approximate_pattern_match(kmer, genome, d).split())


def get_kmer_counts(genome, k, d=0, rc=False):
    """
    Find counts for all approximate kmers
    :param genome: genome string
    :param k: k-mer size
    :param d: maximum distance value
    :param rc: toggle between whether reverse compliments count toward
        the final count
    :return: {approximate_kmer: count}
    """
    approximate_kmers = {}

    # Set initial kmer counts
    kmers = get_all_kmers(k)
    i = 0
    while i <= len(genome) - int(k):
        kmers[genome[i:i + int(k)]] = kmers.get(genome[i:i + int(k)], 0) + 1
        i += 1

    # Handle adding approximate counts if applicable
    for kmer in list(kmers.keys()):
        approximate_kmers[kmer] = approximate_kmers.get(kmer, 0) + kmers[kmer]  # Add count of kmer itself
        for other_kmer in [k for k in kmers.keys() if k != kmer]:
            # For efficiency only run hamming distances if either kmer has a count
            if kmers[kmer] or kmers[other_kmer]:
                # If hamming distance match, add the count of each to the others final count
                if hamming_distance(kmer, other_kmer) <= d:
                    approximate_kmers[kmer] = approximate_kmers.get(kmer, 0) + kmers[other_kmer]
                    approximate_kmers[other_kmer] = approximate_kmers.get(other_kmer, 0) + kmers[kmer]
        # For efficiency in not duplicating hamming distance calculations delete current key
        del kmers[kmer]

    # If reverse compliments is True add counts of the reverse compliment
    if rc:
        check = {}
        for kmer in approximate_kmers:
            # If the kmer is in check, it was already handled by its compliment
            if kmer in check:
                continue
            complement = reverse_complement(kmer)
            approximate_kmers[kmer] = approximate_kmers.get(kmer, 0) + approximate_kmers[complement]
            approximate_kmers[complement] = approximate_kmers[kmer]
            check[complement] = True

    return approximate_kmers


def get_frequent_kmers(genome, k, d=0, rc=False):
    """
    Find the most frequent k-mers for a genome genome.
    :param genome: genome string
    :param k: length of k-mer
    :param d: maximum hamming distance value if approximate method
    :param rc: reverse compliment processing
    :param approx: whether the kmers are absolute or approximate
    :return: space delimited list of the most frequent kmers
    """
    kmers = get_kmer_counts(genome, int(k), int(d), rc)
    value = kmers[max(kmers, key=kmers.get)]
    return ' '.join(sorted([k for k, v in kmers.items() if v == value]))


if __name__ == "__main__":
    # lines = sys.stdin.read().splitlines()
    print(count_approximate_pattern_matches('CCC', 'CATGCCATTCGCATTGTCCCAGTGA', 2))
