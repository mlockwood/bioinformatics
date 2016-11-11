__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


import sys


def pattern_count(text, pattern):
    """
    Count the number of times a pattern appears in a genome text.
    :param text: genome string
    :param pattern: substring value to be searched for within the genome
    :return: count of pattern in genome text
    """
    count = 0
    i = 0
    while i <= len(text) - len(pattern):
        if text[i:i+len(pattern)] == pattern:
            count += 1
        i += 1
    return count


def frequent_kmers(text, k):
    """
    Find the most frequent k-mers for a genome text. A k-mer is a
    string of nucleotides of `k` length. This function will inventory
    all k-mer values within a text and return the values found most
    :param text: genome text
    :param k: length of nucleotide string
    :return: most frequent k-mer nucleotide strings
    """
    frequent_patterns = {}
    i = 0
    while i <= len(text) - int(k):
        frequent_patterns[text[i:i+int(k)]] = frequent_patterns.get(text[i:i+int(k)], 0) + 1
        i += 1
    value = frequent_patterns[max(frequent_patterns, key=frequent_patterns.get)]
    return ' '.join(sorted([k for k, v in frequent_patterns.items() if v == value]))


def pattern_match(substring, text):
    """
    Return a series of indices in which a k-mer substring appears.
    :param substring: string of nucleotides that is being searched for
    :param text: genome text
    :return: all indices of where the substring is found
    """
    index = []
    i = 0
    while i <= len(text) - len(substring):
        if text[i:i+len(substring)] == substring:
            index += [i]
        i += 1
    return ' '.join(index)


def clump_finding(genome, k, L, t):
    """
    Identify all k-mer nucleotide strings that are found of at least
    frequency `t` within a window size (string distance) of `L`.
    :param genome: genome text
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


#lines = sys.stdin.read().splitlines()
#print(clump_finding(*lines, 9, 500, 3)))
