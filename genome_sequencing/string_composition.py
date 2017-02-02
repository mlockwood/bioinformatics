#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

from generic import *


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def string_composition(k, text):
    """
    Produce a set of all k-mer nucleotide sequences in the text.
    :param k: length of nucleotide strings.
    :param text: genome text
    :return: a set of all k-mer nucleotides for the text
    """
    kmers = {}
    i = 0
    k = int(k)
    while i <= (len(text) - k):
        kmers[text[i:i+k]] = True
        i += 1
    return '\n'.join(sorted(list(kmers.keys())))


def string_reconstruction_pre_ordered(kmers):
    """
    Construct a genome path from previously ordered kmers.
    :param kmers: space delimited and ordered kmers for a genome
    :return: the genome path for the kmers
    """
    kmers = kmers.rstrip().split()
    genome_path = kmers[0]
    for kmer in kmers[1:]:
        genome_path += kmer[-1]
    return genome_path


def get_prefix(kmer):
    """
    Return the prefix of a kmer.
    :param kmer: string of nucleotides
    :return: prefix of the kmer
    """
    return kmer[:-1]


def get_suffix(kmer):
    """
    Return the suffix of a kmer.
    :param kmer: string of nucleotides
    :return: suffix of the kmer
    """
    return kmer[1:]


def get_overlap_graph(kmers):
    """
    Construct a genome graph where kmer -> next kmer.
    :param kmers: space delimited and ordered kmers for a genome
    :return: the genome graph for the kmers
    """
    kmers = kmers.rstrip().split()
    graph = {}
    for i in kmers:
        graph[i] = {}
        for j in kmers:
            if get_suffix(i) == get_prefix(j):
                graph[i][j] = True

    return print_overlap_graph(graph)


def print_overlap_graph(graph):
    """
    Print overlap graph by ordered keys where i -> j.
    :param graph: an overlap graph
    :return: printable version of the graph in list representation
    """
    result = ''
    for i in sorted(list(graph.keys())):
        for j in sorted(list(graph[i].keys())):
            result += '{} -> {}\n'.format(i, j)
    return result


def de_bruijn_graph_from_string(k, text):
    """
    Construct a de Bruijn graph from a genome text.
    :param k: the k-mer size
    :param text: genome text
    :return: the de Bruijn graph by list representation
    """
    graph = {}
    node_size = int(k) - 1
    prev = text[0:node_size]
    i = 1
    while i <= (len(text) - node_size):
        if prev not in graph:
            graph[prev] = []
        # Append to node onto from list based on next k-mer sequence
        graph[prev] = graph.get(prev) + [text[i:i+node_size]]
        prev = text[i:i+node_size]
        i += 1

    return print_de_bruijn_graph(graph)


def de_bruijn_graph_by_composition(kmers):
    """
    Construct a de Bruijn graph by composition from individual k-mers.
    :param kmers: space delimited and ordered kmers for a genome
    :return: the de Bruijn graph for the kmers
    """
    kmers = kmers.rstrip().split()
    graph = {}
    for kmer in kmers:
        if get_prefix(kmer) not in graph:
            graph[get_prefix(kmer)] = []
        graph[get_prefix(kmer)] = graph.get(get_prefix(kmer)) + [get_suffix(kmer)]

    return print_de_bruijn_graph(graph)


def print_de_bruijn_graph(graph):
    """
    Print de Bruijn graph by ordered keys where i -> j.
    :param graph: a de Bruijn graph
    :return: printable version of the graph in list representation
    """
    result = ''
    for i in sorted(list(graph.keys())):
        result += '{} -> {}\n'.format(i, ','.join(sorted(graph[i])))
    return result


def eulerian_path(graph, nearly=False):
    """
    Traverse a Eulerian graph to find a cycle where each edge is
    visited exactly once.
    :param graph: {node: {to_nodes: True}}
    :param nearly: if graph is nearly balanced will find Eulerian path
    :return: string representation of a path
    """
    if nearly:
        # Count the balance of edges at each node
        counts = {}
        for node in graph:
            # Add for the from node
            if node not in counts:
                counts[node] = {'from': len(graph[node]), 'to': 0}
            else:
                counts[node]['from'] = len(graph[node])

            # Add for each to node
            for to_node in graph[node]:
                if to_node not in counts:
                    counts[to_node] = {'from': 0, 'to': 1}
                else:
                    counts[to_node]['to'] = counts[to_node].get('to', 0) + 1

        # Find the missing edge, note it then add to graph
        for node in counts:
            if counts[node]['from'] < counts[node]['to']:
                from_node = node
            elif counts[node]['from'] > counts[node]['to']:
                to_node = node
        final_edge = [from_node, to_node]
        if from_node not in graph:
            graph[from_node] = {}
        graph[from_node][to_node] = True

    cycle = [next(iter(graph))]

    while graph:
        # If there are more to_nodes for the last node in the cycle, continue
        if cycle[-1] in graph:
            cycle.append(next(iter(graph[cycle[-1]])))
            del graph[cycle[-2]][cycle[-1]]

            # Remove nodes completely if they have no more to_nodes
            if not graph[cycle[-2]]:
                del graph[cycle[-2]]

        # If there are no more to_nodes choose a new_start location
        else:
            i = 0
            while i < len(cycle):
                # If a node has to_nodes restructure cycle to start at this node
                if cycle[i] in graph:
                    cycle = cycle[i:-1] + cycle[:i+1]
                    i = len(cycle)
                i += 1

            # If no new start was chosen the graph is not Eulerian
            if i == len(cycle):
                raise ValueError('Graph is not Eulerian.')

    # Handle realignment of nearly balanced Eulerian path
    if nearly:
        i = 0
        while i < len(cycle):
            if cycle[i:i+2] == final_edge:
                cycle = cycle[i+1:-1] + cycle[:i+1]
                i = len(cycle)
            i += 1

    return '->'.join(cycle)


lines = sys.stdin.read().splitlines()
print(eulerian_path(lines_to_graph_dict(lines), nearly=True))
