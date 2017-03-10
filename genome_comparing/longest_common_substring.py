#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

from constants import *
from generic import *


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def solve_LCS(v, w, method='default', sigma=0, mu=0):
    """
    Run and output results for the LCS problem.
    :param v: first string
    :param w: second string
    :param method: the method for scoring
    :param sigma: penalty for indels
    :param mu: match and mismatch scoring
    :return: a LCS backtrack tracer
    """
    # For fitting run all substrings of `v` against `w`
    if re.search('fitting', method):
        best_score = (0, None, '')  # (score, backtrack, v_substring)
        i = 0
        while i < len(v):
            score, backtrack, pointer = run_LCS(v[i:], w, method, sigma, mu)
            if score > best_score[0]:
                best_score = (score, backtrack, v[i:i+pointer[0]+1])
            i += 1

    # Otherwise perform a standard run
    else:
        score, backtrack, pointer = run_LCS(v, w, method, sigma, mu)

    # Prepare outputs
    if method == 'global':
        return score, output_LCS_alignment(backtrack, v, w, len(v), len(w))

    elif method == 'local':
        return score, output_LCS_alignment(backtrack, v, w, pointer[0], pointer[1])

    elif re.search('fitting', method):
        return best_score[0], output_LCS_alignment(best_score[1], best_score[2], w, len(best_score[2]), len(w))

    elif method == 'edit_distance':
        return -score

    else:
        return output_LCS(backtrack, v, len(v), len(w))


def run_LCS(v, w, method, sigma, mu):
    """
    Create a LCS backtrack tracer for strings v and w.
    :param v: first string
    :param w: second string
    :param method: the method for scoring
    :param sigma: penalty for indels
    :param mu: match and mismatch scoring
    :return: a LCS backtrack tracer
    """
    path = {0: {0: 0}}
    multiplier = -sigma if method != 'edit_distance' else -1

    # Initialize the first column
    i = 1
    while i <= len(v):
        path[i] = {0: i * multiplier}
        i += 1

    # Initialize the first row
    j = 1
    while j <= len(w):
        path[0][j] = j * multiplier
        j += 1

    # Process each subsequent row
    backtrack = {}
    best_score = (-sys.maxsize, None)
    i = 1
    while i <= len(v):
        backtrack[i] = {}
        j = 1
        while j <= len(w):
            # Build the options list depending on the method selected
            if method == 'global' or method == 'fitting-global':
                options = [
                    (path[i-1][j] - sigma, 'D'),
                    (path[i][j-1] - sigma, 'I'),
                    (path[i-1][j-1] + mu[v[i-1]][w[j-1]], 'M')
                ]

            elif method == 'local':
                options = [
                    (0, 'F'),
                    (path[i-1][j] - sigma, 'D'),
                    (path[i][j-1] - sigma, 'I'),
                    (path[i-1][j-1] + mu[v[i-1]][w[j-1]], 'M')
                ]

            elif method == 'edit_distance':
                options = [
                    (path[i-1][j] - 1, 'D'),
                    (path[i][j-1] - 1, 'I'),
                    (path[i-1][j-1], 'M') if v[i-1] == w[j-1] else (path[i-1][j-1] - 1, 'U')
                ]

            else:
                options = [
                    (path[i-1][j] - sigma, 'D'),
                    (path[i][j-1] - sigma, 'I'),
                    (path[i-1][j-1] + 1, 'M') if v[i-1] == w[j-1] else (path[i-1][j-1] - mu, 'U')
                ]

            # Set path and backtrack based on which option has the maximum score
            path[i][j], backtrack[i][j] = max(options)

            # Handle best scores
            if path[i][j] > best_score[0] and (method == 'local' or method == 'fitting'):
                best_score = (path[i][j], (i, j))
            j += 1
        i += 1

    return best_score[0] if best_score[0] != -sys.maxsize else path[i-1][j-1], backtrack, best_score[-1]


def output_LCS(backtrack, v, i, j):
    """
    Using a LCS backtrack and the input string construct the LCS.
    :param backtrack: LCS backtrack tracer
    :param v: first string
    :param i: length of v
    :param j: length of w, the second string
    :return: LCS
    """
    if i == 0 or j == 0:
        return ''
    if backtrack[i][j] == 'D':
        return output_LCS(backtrack, v, i-1, j)
    elif backtrack[i][j] == 'I':
        return output_LCS(backtrack, v, i, j-1)
    elif backtrack[i][j] == 'U':
        return output_LCS(backtrack, v, i-1, j-1)
    elif backtrack[i][j] == 'M':
        return output_LCS(backtrack, v, i-1, j-1) + v[i-1]


def output_LCS_alignment(backtrack, v, w, i, j):
    """
    Using a LCS backtrack and both input strings construct an alignment
    between the two string.
    :param backtrack: LCS backtrack tracer
    :param v: first string
    :param w: second string
    :param i: length of v
    :param j: length of w, the other string
    :return: LCS
    """
    # Base case, consider if there are any remaining elements from v or w
    if i == 0 or j == 0:
        if i > 0:
            return v[:i], '-' * i
        elif j > 0:
            return '-' * j, w[:j]
        else:
            return '', ''

    # Base care for local - free ride to the node
    if backtrack[i][j] == 'F':
        return '', ''

    # Deletion
    if backtrack[i][j] == 'D':
        new_v, new_w = output_LCS_alignment(backtrack, v, w, i-1, j)
        return new_v + v[i-1], new_w + '-'

    # Insertion
    elif backtrack[i][j] == 'I':
        new_v, new_w = output_LCS_alignment(backtrack, v, w, i, j-1)
        return new_v + '-', new_w + w[j-1]

    # Match or mismatch
    elif backtrack[i][j] == 'M' or backtrack[i][j] == 'U':
        new_v, new_w = output_LCS_alignment(backtrack, v, w, i-1, j-1)
        return new_v + v[i-1], new_w + w[j-1]


class DAG(object):

    def __init__(self, graph, root, tail):
        self.graph = graph
        self.root = root
        self.tail = tail
        self.longest_length = 0
        self.longest_path = []

    @staticmethod
    def build_graph(lines):
        return lines_to_graph_dict(lines, weighted=True)

    def longest_path_in_DAG(self):
        """
        Take a directed acyclic graph (DAG) and find the longest path by
        weight.
        :return: length and path
        """
        if not self.longest_path:
            self.explore_DAG([self.root])
        return self.longest_length, self.longest_path

    def explore_DAG(self, path, weight=0):
        """
        Explore a DAG recursively seeking for the best solution.
        :param path: previously explored path
        :param weight: current weight given path
        :return: save results to self.longest_[length|path]
        """

        # Handle base case where the tail has been reached
        if path[-1] == self.tail:
            if weight > self.longest_length:
                self.longest_length = weight
                self.longest_path = path

        # Otherwise explore other nodes
        if path[-1] in self.graph:
            current_tail = path[-1]
            for node in self.graph[path[-1]]:
                self.explore_DAG(path+[node], weight+self.graph[current_tail][node])


lines = sys.stdin.read().splitlines()
score, out = solve_LCS(*lines, method='fitting', sigma=1, mu=1)
print(score)
print(out[0])
print(out[1])
