#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import copy
import sys

from constants import *
from generic import *


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def solve_LCS(v, w, method='default', sigma=0, mu=0, epsilon=0):
    """
    Run and output results for the LCS problem.
    :param v: first string
    :param w: second string
    :param method: the method for scoring
    :param sigma: penalty for indels (particularly those starting gaps)
    :param mu: match and mismatch scoring
    :param epsilon: penalty for indels after a gap has been started
    :return: a LCS backtrack tracer
    """
    path, backtrack, score, i, j = run_LCS(v, w, method, sigma, mu, epsilon)

    # Prepare outputs
    if method == 'global':
        return path[len(v)][len(w)], output_LCS_alignment(backtrack, v, w, len(v), len(w))

    elif method == 'local':
        return path[i][j], output_LCS_alignment(backtrack, v, w, i, j)

    elif re.search('fitting', method):
        return path[i][j], output_LCS_alignment(backtrack, v[:i], w, i, j, method)

    elif re.search('overlap', method):
        return path[i][j], output_LCS_alignment(backtrack, v, w[:j], i, j, method)

    elif method == 'affine-gaps':
        return path[len(v)][len(w)], output_LCS_alignment(backtrack, v, w, len(v), len(w))

    elif method == 'edit_distance':
        return -path[len(v)][len(w)]

    else:
        return output_LCS(backtrack, v, len(v), len(w))


def run_LCS(v, w, method, sigma, mu, epsilon):
    """
    Create a LCS backtrack tracer for strings v and w.
    :param v: first string
    :param w: second string
    :param method: the method for scoring
    :param sigma: penalty for indels (particularly those starting gaps)
    :param mu: match and mismatch scoring
    :param epsilon: penalty for indels after a gap has been started
    :return: a LCS backtrack tracer
    """
    path = {0: {0: 0}}
    if method == 'edit_distance':
        multiplier = -1
    elif re.search('fitting', method):
        multiplier = 0
    else:
        multiplier = -sigma

    # Initialize the first column
    i = 1
    while i <= len(v):
        path[i] = {0: i * multiplier if not re.search('overlap', method) else 0}
        i += 1

    # Initialize the first row
    j = 1
    while j <= len(w):
        path[0][j] = j * multiplier
        j += 1

    lower, upper = copy.deepcopy(path), copy.deepcopy(path)

    # Process each subsequent row
    backtrack = {}
    best_score = (-sys.maxsize, None, None)
    i = 1
    while i <= len(v):
        backtrack[i] = {}
        j = 1
        while j <= len(w):
            # Handle mu value
            if isinstance(mu, dict):
                mu_value = mu[v[i-1]][w[j-1]]
            else:
                mu_value = -mu if v[i-1] != w[j-1] else 1

            # Build the options list depending on the method selected
            if method == 'global' or method == 'fitting-global':
                options = [
                    (path[i-1][j] - sigma, 'D'),
                    (path[i][j-1] - sigma, 'I'),
                    (path[i-1][j-1] + mu_value, 'M')
                ]

            elif method == 'local':
                options = [
                    (0, 'F'),
                    (path[i-1][j] - sigma, 'D'),
                    (path[i][j-1] - sigma, 'I'),
                    (path[i-1][j-1] + mu_value, 'M')
                ]

            elif method == 'affine-gaps':
                # Set lower scores
                lower[i][j] = max([
                    lower[i-1][j] - epsilon,
                    path[i-1][j] - sigma
                ])

                # Set upper scores
                upper[i][j] = max([
                    upper[i][j-1] - epsilon,
                    path[i][j-1] - sigma
                ])

                # Set options for path ("middle")
                options = [
                    (lower[i][j], 'D'),
                    (upper[i][j], 'I'),
                    (path[i-1][j-1] + mu_value, 'M')
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
                    (path[i-1][j-1] + mu_value, 'M' if v[i-1] == w[j-1] else 'U')
                ]

            # Set path and backtrack based on which option has the maximum score
            path[i][j], backtrack[i][j] = max(options)

            # Handle best scores
            if path[i][j] > best_score[0]:
                if (method == 'local' or
                        (j == len(w) and re.search('fitting', method)) or
                        (i == len(v) and re.search('overlap', method))):
                    best_score = (path[i][j], i, j)
            j += 1
        i += 1

    return path, backtrack, best_score[0], best_score[1], best_score[2]


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


def output_LCS_alignment(backtrack, v, w, i, j, method=''):
    """
    Using a LCS backtrack and both input strings construct an alignment
    between the two string.
    :param backtrack: LCS backtrack tracer
    :param v: first string
    :param w: second string
    :param i: length of v
    :param j: length of w, the other string
    :param method: method check for allowing extra indels at BOS/EOS
    :return: LCS
    """
    # Base case, consider if there are any remaining elements from v or w
    if i == 0 or j == 0:
        if i > 0 and not re.search('fitting|overlap', method):
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
        new_v, new_w = output_LCS_alignment(backtrack, v, w, i-1, j, method)
        return new_v + v[i-1], new_w + '-'

    # Insertion
    elif backtrack[i][j] == 'I':
        new_v, new_w = output_LCS_alignment(backtrack, v, w, i, j-1, method)
        return new_v + '-', new_w + w[j-1]

    # Match or mismatch
    elif backtrack[i][j] == 'M' or backtrack[i][j] == 'U':
        new_v, new_w = output_LCS_alignment(backtrack, v, w, i-1, j-1, method)
        return new_v + v[i-1], new_w + w[j-1]


def score_alignments(v, w, match=1, sigma=1, mu=1):
    """
    Run and output results for the LCS problem.
    :param v: first string
    :param w: second string
    :param match: the points awarded for a match
    :param sigma: penalty for indels
    :param mu: the points lost for a mismath
    :return: score
    """
    score = 0
    i = 0
    while i < len(v):
        # If a match, add match to score
        if v[i] == w[i]:
            score += match

        # If an indel has been found (for insert or delete) subtract sigma from score
        elif v[i] == '-' or w[i] == '-':
            score -= sigma

        # If a mismatch, subtract mu from score
        else:
            score -= mu
        i += 1
    return score


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
print(solve_LCS(*lines, method='affine-gaps', sigma=11, mu=BLOSUM62, epsilon=1))
