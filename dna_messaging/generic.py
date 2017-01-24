#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dna_messaging.constants import BASES


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def get_all_kmers(k, kmer=''):
    """
    Generate a dictionary of all kmers for k
    :param k: length of kmer
    :param kmer: a dummy value that will be built into kmers
    :return: {kmer: 0}
    """
    kmers = {}
    # Base case by adding the kmer to resulting output
    if len(kmer) == k:
        return {kmer: 0}

    # Recurse another layer (length) for each base in BASES
    else:
        for base in BASES:
            kmers.update(get_all_kmers(k, '{}{}'.format(kmer, base)))
        return kmers


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


def get_all_mismatched_kmers(pattern, d, kmer=''):
    """
    Generate a dictionary of all mismatched kmers for the pattern who
    have a hamming distance less than or equal to d.
    :param pattern: mismatched kmers will be generated from this
    :param d: maximum hamming distance between a kmer and the pattern
    :param kmer: a dummy value that will be built into kmers
    :return: {kmer: hamming distance}
    """
    kmers = {}
    # Base case by adding the kmer to resulting output
    if len(kmer) == len(pattern):
        return {kmer: hamming_distance(pattern, kmer)}

    # Recurse another layer (length) for each base in BASES
    else:
        # If the maximum amount of distance has been reached, only pursue the base that matches pattern for index
        if hamming_distance(pattern[:len(kmer)], kmer) == int(d):
            kmers.update(get_all_mismatched_kmers(pattern, d, '{}{}'.format(kmer, pattern[len(kmer)])))

        # Otherwise seek all base options
        else:
            for base in BASES:
                kmers.update(get_all_mismatched_kmers(pattern, d, '{}{}'.format(kmer, base)))
        return kmers
