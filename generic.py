#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re

from constants import BASES


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def reverse_complement(string):
    lookup = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    text = ''
    i = -1
    while i >= -len(string):
        text += lookup[string[i]]
        i -= 1
    return text


def get_all_binary_kmers(k, kmer=''):
    """
        Generate a dictionary of all kmers for k. This is for the
        universal string problem so these are binary.
        :param k: length of kmer
        :param kmer: a dummy value that will be built into kmers
        :return: {kmer: 0}
        """
    kmers = {}
    # Base case by adding the kmer to resulting output
    if len(kmer) == k:
        return {kmer: 0}

    # Recurse another layer (length) for each base in [0, 1]
    else:
        for base in ['0', '1']:
            kmers.update(get_all_binary_kmers(k, '{}{}'.format(kmer, base)))
        return kmers


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


def lines_to_graph_dict(lines):
    """
    Take a graph composed of line delimited x -> y, z node to node(s)
    and convert this to a dictionary representation of
    {node: {to_node: count}}.
    :param lines: list of string representations of the graph nodes
    :return: graph dictionary
    """
    graph = {}
    for line in lines:
        # Split the line into the first node left of the -> and a list of to_nodes from the right
        node, to_nodes = re.split('->', re.sub('\s', '', line))

        # If the node has not been seen already add it to the graph with a value of an empty dict
        if node not in graph:
            graph[node] = {}

        # Set each to_node for the node
        for to_node in re.split(',', to_nodes):
            graph[node][to_node] = graph[node].get(to_node, 0) + 1

    return graph


def select_scores_with_ties(scores, n):
    """
    Return the N best elements in a score after factoring in tie
    breaking. Scores should be sorted before function use.
    :param scores: [(score, element), (score, element)]
    :return: new N index
    """
    # Change N in the event of a tie
    i = int(n)
    while i < len(scores):
        if scores[i][0] == scores[i - 1][0]:
            i += 1
        else:
            break
    return i