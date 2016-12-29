__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


import sys


def LCS_back_track(v, w):
    """
    Create a LCS back track tracer for strings v and w.
    :param v: first string
    :param w: second string
    :return: a LCS back track tracer
    """
    path = {0: {0: 0}}
    backtrack = {}

    # Initialize the first column
    i = 0
    while i <= len(v):
        path[i] = {0: 0}
        i += 1

    # Initialize the first row
    j = 0
    while j <= len(w):
        path[0][j] = 0
        j += 1

    # Process each subsequent row
    i = 1
    while i <= len(v):
        backtrack[i] = {}
        j = 1
        while j <= len(w):
            diag = path[i-1][j-1] + 1 if v[i-1] == w[j-1] else path[i-1][j-1]
            path[i][j] = max([path[i-1][j], path[i][j-1], diag])

            # Set pointers based on which option had the maximum value from above
            if path[i][j] == path[i-1][j]:
                backtrack[i][j] = 'G'
            elif path[i][j] == path[i][j-1]:
                backtrack[i][j] = 'B'
            elif path[i][j] == (path[i-1][j-1] + 1) and v[i-1] == w[j-1]:
                backtrack[i][j] = 'R'
            j += 1
        i +=1

    return backtrack


def output_LCS(backtrack, v, i, j):
    if i == 0 or j == 0:
        return ''
    if backtrack[i][j] == 'G':
        return output_LCS(backtrack, v, i-1, j)
    elif backtrack[i][j] == 'B':
        return output_LCS(backtrack, v, i, j-1)
    elif backtrack[i][j] == 'R':
        return output_LCS(backtrack, v, i-1, j-1) + v[i-1]


def solve_LCS(v, w):
    backtrack = LCS_back_track(v, w)
    return output_LCS(backtrack, v, len(v), len(w))


sys.setrecursionlimit(10000)
print(solve_LCS(*sys.stdin.read().splitlines()))
