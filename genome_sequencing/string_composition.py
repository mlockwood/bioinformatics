#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import copy
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
    i = 1
    while i <= (len(text) - node_size):
        if text[i-1:i-1+node_size] not in graph:
            graph[text[i-1:i-1+node_size]] = {}
        graph[text[i-1:i-1+node_size]][text[i:i+node_size]] = graph[text[i-1:i-1+node_size]].get(text[i:i+node_size],
                                                                                                 0) + 1
        i += 1

    return print_de_bruijn_graph(graph)


def de_bruijn_graph_by_composition(kmers):
    """
    Construct a de Bruijn graph by composition from individual k-mers.
    :param kmers: list of k-mers for a genome
    :return: the de Bruijn graph for the kmers
    """
    graph = {}
    for kmer in kmers:
        if get_prefix(kmer) not in graph:
            graph[get_prefix(kmer)] = {}
        graph[get_prefix(kmer)][get_suffix(kmer)] = graph[get_prefix(kmer)].get(get_suffix(kmer), 0) + 1

    return graph


def print_de_bruijn_graph(graph):
    """
    Print de Bruijn graph by ordered keys where i -> j.
    :param graph: a de Bruijn graph
    :return: printable version of the graph in list representation
    """
    result = ''
    for i in sorted(list(graph.keys())):
        to_nodes = []
        for key in sorted(graph[i].keys()):
            for x in range(graph[i][key]):
                to_nodes.append(key)
        result += '{} -> {}\n'.format(i, ','.join(to_nodes))
    return result


def count_graph_edges(graph):
    """
    Traverse a graph where {from_node: {to_node: count}} and produce a
    dictionary where {node: {from: X, to: X} where X is the count of
    how many edges come from or to the node.
    :param graph: a graph of nodes to nodes
    :return: counts of from and to edges for each node
    """
    # Count the balance of edges at each node
    counts = {}
    lookup = {}
    for node in graph:
        for to_node in graph[node]:
            # Handle from counts first
            if node not in counts:
                counts[node] = {'from': graph[node][to_node], 'to': 0}
            else:
                counts[node]['from'] = counts[node].get('from', 0) + graph[node][to_node]

            # Then handle to counts
            if to_node not in counts:
                counts[to_node] = {'from': 0, 'to': graph[node][to_node]}
            else:
                counts[to_node]['to'] = counts[to_node].get('to', 0) + graph[node][to_node]

            # Build a lookup in the event that the to_node has exactly one inbound and one outbound edge
            lookup[to_node] = node

    return counts, lookup


def set_guaranteed_paths(graph, counts, lookup):
    """
    Find all guaranteed (maximal non-branching) paths and then remove
    them from the graph attaching them to their root/from node.
    :param graph: a graph of nodes to nodes
    :param counts: counts of from and to edges for each node
    :param lookup: reference for connecting a to node to its from node
    :return: modified graph with guaranteed paths and root if it exists
    """
    # Determine if there is a definitive root instance and if there are any guaranteed paths
    root = None
    guaranteed = {}
    for node in counts:
        if counts[node]['from'] == 1 and counts[node]['to'] == 1:
            guaranteed[node] = next(iter(graph[node]))
        if counts[node]['from'] > counts[node]['to']:
            root = node

    # Remove guaranteed paths from graph
    for node in guaranteed:
        del graph[node]

    # Transfer guaranteed paths to ordered lists within the graph DS
    final = {}
    prev = []
    key = None
    prev_key = None
    while guaranteed:
        # Handle key assignments
        if key == prev_key:
            key = next(iter(guaranteed))

        # Handle building the path
        if not prev:
            prev.append(key)
        prev.append(guaranteed[key])
        prev_key = key

        # First test if the value is in guaranteed (ergo it has not been selected yet)
        if guaranteed[key] in guaranteed:
            key = guaranteed[key]

        # Next test if the value is in the final (ergo it has already been selected)
        elif guaranteed[key] in final:
            final[prev[0]] = prev[1:] + final[guaranteed[key]]
            del final[guaranteed[key]]
            prev = []

        # Add to final if the path has finished and will not continue guaranteed
        else:
            final[prev[0]] = prev[1:]
            prev = []

        del guaranteed[prev_key]

    # To prepare the final graph, take the final DS of guaranteed paths and use the lookup to move these to the graph
    for node in final:
        if lookup[node] in graph:
            graph[lookup[node]][node] = final[node]
        elif node == final[node][-1]:
            if node not in graph:
                graph[node] = {}
            graph[node][final[node][0]] = final[node][1:]

    return graph, root


def eulerian_path(graph, nearly=False):
    """
    Traverse a Eulerian graph to find a cycle where each edge is
    visited exactly once.
    :param graph: {node: {to_nodes: True}}
    :param nearly: if graph is nearly balanced will find Eulerian path
    :return: string representation of a path
    """
    if nearly:
        counts = count_graph_edges(graph)[0]

        # Find the missing edge, note it then add to graph
        from_node, to_node = (None, None)
        for node in counts:
            if counts[node]['from'] < counts[node]['to']:
                from_node = node
            elif counts[node]['from'] > counts[node]['to']:
                to_node = node

        # Test to be sure that the graph is nearly balanced
        if not from_node and not to_node:
            nearly = False
        else:
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


def resolve_overlaps(graph_string):
    """
    Take overlapping k-mers in a graph_string and combine them
    :param graph_string: result of eulerian_path
    :return: string composition by removing overlap
    """
    kmers = re.split('->', graph_string)
    composition = kmers[0]
    for kmer in kmers[1:]:
        composition += kmer[-1]
    return composition


def genome_reconstruction(kmers):
    """
    Also known as the string reconstruction problem, take k-mers and
    reassemble the string they are supposed to represent.
    :param kmers: list of k-mers
    :return: genome reconstructed from the k-mers
    """
    de_bruijn = de_bruijn_graph_by_composition(kmers)
    graph_string = eulerian_path(de_bruijn, nearly=True)
    return resolve_overlaps(graph_string)


def universal_circular_string(k):
    """
    Find a universal circular string for binary k-mers of length k.
    :param k: length of k-mers
    :return: universal string
    """
    return genome_reconstruction(list(get_all_binary_kmers(int(k)).keys()))[:-int(k)+1]


def strings_spelled_by_gapped_patterns(pairs, k, d):
    """
    Take read pairs for a genome of length k and read distance d and
    test whether the genome string for them is viable.
    :param pairs: read pairs
    :param k: length of k-mers
    :param d: read distance
    :return: string spelled by read pairs or no finding result
    """
    # Build the prefix and suffix strings from the read pairs
    prefix, suffix = re.split('\|', pairs[0])
    for pair in pairs[1:]:
        pre, suf = re.split('\|', pair)
        prefix += pre[-1]
        suffix += suf[-1]

    # Compare each position
    i = int(k) + int(d)
    while i < len(prefix):
        if prefix[i] != suffix[i-int(k)-int(d)]:
            return "there is no string spelled by the gapped patterns"
        i += 1
    return prefix + suffix[-int(k)-int(d):]


def read_pairs_path(pairs, k, d):
    """
    Take read pairs for a genome of length k and read distance d
    and find their correct ordering and resulting genome string.
    :param pairs: read pairs
    :param k: length of k-mers
    :param d: read distance
    :return: genome string that fits the read pairs
    """
    # Build the de Bruijn graph for the read pairs
    graph = read_pairs_de_bruijn_graph_by_composition(pairs)
    counts, lookup = count_graph_edges(graph)
    graph, root = set_guaranteed_paths(graph, counts, lookup)

    # Convert remaining graph items for search efficiency
    for node in graph:
        new_dict = {}
        for to_node in graph[node]:
            if to_node[0] not in new_dict:
                new_dict[to_node[0]] = {}
            new_dict[to_node[0]][to_node[1]] = graph[node][to_node]
        graph[node] = copy.deepcopy(new_dict)

    # If there is a root only search paths that start at the root, otherwise allow any node to be the root
    if root:
        return read_pair_genome_reconstruction(search_paired_de_bruijn([root], graph, int(k), int(d)), k, d)
    else:
        for node in graph:
            new_path = search_paired_de_bruijn([node], graph, int(k), int(d))
            if new_path:
                return read_pair_genome_reconstruction(new_path, k, d)


def read_pairs_de_bruijn_graph_by_composition(pairs):
    """
    Take a collection of kmer read pairs delimited by | and construct a
    de Bruijn graph from them.
    :param pairs: read pairs of k-mers
    :return: de Bruijn graph in dictionary form
    """
    graph = {}
    for kmer_tuple in pairs:

        # Split the string into a list of kmers
        kmer_tuple = re.split('\|', kmer_tuple)

        # Create tuples of the prefixes and suffixes of the kmers
        prefix = tuple([get_prefix(kmer) for kmer in kmer_tuple])
        suffix = tuple([get_suffix(kmer) for kmer in kmer_tuple])

        # Set the graph to {prefix: {suffix: count}}
        if prefix not in graph:
            graph[prefix] = {}
        graph[prefix][suffix] = graph[prefix].get(suffix, 0) + 1

    return graph


def search_paired_de_bruijn(path, graph, k, d):
    """
    This function manages the actual searching for appropriate paths.
    It has a few built-in strategies to cut down run time; allowing a
    root start, requiring the next path's start pair is equal to the
    final pair of index -k-d in the previous path, and finding
    guaranteed paths (nodes with exactly one inbound and one outbound
    edge) which significantly reduces recursion depth.
    :param path: previous path
    :param graph: remaining unexplored edges
    :param k: length of k-mers
    :param d: read distance
    :return: root path if graph has been fully explored
    """
    # Base case is an empty graph
    if not graph:
        return path

    # First test that there are still next options given the last explored node
    if path[-1] in graph:
        # If existing path has at least k+d nodes, ensure that the next node matches the second k-mer of the pair at
        # path index -k-d
        if len(path) > k + d:
            if path[-k-d][1] in graph[path[-1]]:
                for to_node1 in graph[path[-1]][path[-k-d][1]]:
                    out = delete_used_node_dependencies(graph, path, (path[-k-d][1], to_node1), k, d)
                    if out:
                        return out
        else:
            print('\n', path[-1], graph[path[-1]])
            for to_node0 in graph[path[-1]]:
                print(path)
                for to_node1 in graph[path[-1]][to_node0]:
                    out = delete_used_node_dependencies(graph, path, (to_node0, to_node1), k, d)
                    if out:
                        return out


def delete_used_node_dependencies(graph, path, to_node, k, d):
    """
    This creates a copy of the graph and deletes/removes the node that
    was just used.
    :param path: previous path
    :param graph: remaining unexplored edges + recent most explored
    :param to_node: currently selected node -> (k-mer, k-mer)
    :param k: length of k-mers
    :param d: read distance
    :return: graph copy less any dependencies
    """
    next_graph = copy.deepcopy(graph)
    from_node = path[-1]
    path.append(to_node)

    # Check if the next node is a guaranteed path
    if isinstance(next_graph[from_node][to_node[0]][to_node[1]], list):
        g_path = next_graph[from_node][to_node[0]][to_node[1]]
        path += g_path

        # Check the additions in the path to ensure their read pairs match the paradigm
        i = -1
        while i >= -(len(g_path) + 1 if len(g_path) + 1 + k + d < len(path) else len(path) - k - d):
            if path[i][0] != path[i-k-d][1]:
                print('Guaranteed path error with mismatched read pairs')
                return None
            i -= 1

        # This will force deletion of dependencies on next conditional
        next_graph[from_node][to_node[0]][to_node[1]] = 0

    # Test if the count can be reduced
    if next_graph[from_node][to_node[0]][to_node[1]] > 1:
        next_graph[from_node][to_node[0]][to_node[1]] = next_graph[from_node][to_node[0]].get(to_node[1], 0) - 1

    # Otherwise begin deletion of dependencies
    else:
        del next_graph[from_node][to_node[0]][to_node[1]]
        if not next_graph[from_node][to_node[0]]:
            del next_graph[from_node][to_node[0]]
        if not next_graph[from_node]:
            del next_graph[from_node]

    return search_paired_de_bruijn(path, next_graph, k, d)


def read_pair_genome_reconstruction(pairs, k, d):
    """
    Take read pairs for a genome of length k and read distance d and
    construct their genome string. These must already be ordered tuples
    and they must have already been tested to ensure their positions
    match/overlap properly.
    :param pairs: read pairs
    :param k: length of k-mers
    :param d: read distance
    :return: string spelled by read pairs
    """
    # Build the prefix and suffix strings from the read pairs
    if pairs:
        prefix, suffix = pairs[0]
        for pair in pairs[1:]:
            prefix += pair[0][-1]
            suffix += pair[1][-1]

        return prefix + suffix[-int(k)-int(d):]


def contig_generation(kmers):
    """
    From a set of kmers find all contiguous non-branching paths in a de
    Bruijn graph that represent contiguous sections of a genome text.
    :param kmers: list of k-mers for a genome
    :return: all contigs for the k-mers
    """
    paths = maximal_non_branching_paths(de_bruijn_graph_by_composition(kmers))
    final = []
    for path in paths:
        kmer = path[0]
        for node in path[1:]:
            kmer += node[-1]
        final.append(kmer)
    return final


def maximal_non_branching_paths(graph):
    """
    Take a graph of nodes to nodes and find all of the maximal non-
    branching paths.
    :param graph: a graph of nodes to nodes
    :return: all maximal non-branching paths
    """
    counts, lookup = count_graph_edges(graph)
    graph = set_guaranteed_paths(graph, counts, lookup)[0]

    # Take graph and build sorted maximal non-branching paths
    paths = []
    for node in graph:
        for to_node in graph[node]:
            if isinstance(graph[node][to_node], list):
                paths.append((node, to_node, *graph[node][to_node]))
            else:
                for x in range(graph[node][to_node]):
                    paths.append((node, to_node))
    return sorted(paths)
