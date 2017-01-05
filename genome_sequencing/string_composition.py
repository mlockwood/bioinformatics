__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


import sys


def string_composition(k, text):
    """
    Produce a set of all k-mer nucleotide sequences in the text.
    :param k: length of nucleotide strings
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
    Return the prefix of a kmer
    :param kmer: string of nucleotides
    :return: prefix of the kmer
    """
    return kmer[:-1]


def get_suffix(kmer):
    """
    Return the suffix of a kmer
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


#lines = sys.stdin.read()
#print(de_bruijn_graph_by_composition(lines))
