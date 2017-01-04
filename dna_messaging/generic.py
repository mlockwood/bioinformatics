#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dna_messaging.constants import BASES


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


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
