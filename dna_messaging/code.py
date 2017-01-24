import copy
import math
from numpy.random import choice
import random
import sys

BASES = {'A', 'C', 'G', 'T'}
CONVERT = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
PROFILE_INDEX = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}


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
    Generate a dictionary of all mismatched kmers for the pattern that
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


def reverse_complement(string):
    lookup = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    text = ''
    i = -1
    while i >= -len(string):
        text += lookup[string[i]]
        i -= 1
    return text


def pattern_count(genome, pattern):
    """
    Count the number of times a pattern appears in a genome genome.
    :param genome: genome string
    :param pattern: substring value to be searched for within the genome
    :return: count of pattern in genome string
    """
    count = 0
    i = 0
    while i <= len(genome) - len(pattern):
        if genome[i:i+len(pattern)] == pattern:
            count += 1
        i += 1
    return count


def pattern_match(substring, genome):
    """
    Return a series of indices in which a k-mer substring appears.
    :param substring: string of nucleotides that is being searched for
    :param genome: genome string
    :return: all indices of where the substring is found
    """
    index = []
    i = 0
    while i <= len(genome) - len(substring):
        if genome[i:i+len(substring)] == substring:
            index += [i]
        i += 1
    return ' '.join([str(i) for i in index])


def clump_finding(genome, k, L, t):
    """
    Identify all k-mer nucleotide strings that are found of at least
    frequency `t` within a window size (string distance) of `L`.
    :param genome: genome string
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
        if i >= (L-k+1):
            del window[lookup[i-(L-k+1)]][i-(L-k+1)]

        # Add kmer to results if there are at least 't' occurrences within the range
        if len(window[kmer]) >= t:
            results[kmer] = True

        i += 1

    return ' '.join(sorted(list(results.keys())))


def approximate_pattern_match(kmer, genome, d):
    """
    Find the indices where the kmer appears in the genome with a
    hamming distance up to `d`
    :param kmer: kmer string
    :param genome: genome string
    :param d: maximum distance value
    :return: space delimited list of indices
    """
    i = 0
    indices = []
    while i <= len(genome) - len(kmer):
        if hamming_distance(kmer, genome[i:i+len(kmer)]) <= int(d):
            indices.append(i)
        i += 1
    return ' '.join(str(i) for i in indices)


def count_approximate_pattern_matches(kmer, genome, d):
    """
    Find the count of how often the kmer appears in the genome with a
    hamming distance up to `d`
    :param kmer: kmer string
    :param genome: genome string
    :param d: maximum distance value
    :return: count
    """
    return len(approximate_pattern_match(kmer, genome, d).split())


def computing_frequencies(genome, k):
    """
    Frequency array for k-mers in a genome. This is inefficient and is
    not used as a part of any of the algorithms herein.
    :param genome: genome string
    :param k: length of k-mer
    :return: space delimited list of the frequencies for each k-mer
    """
    kmers = get_all_kmers(int(k))
    i = 0
    while i <= len(genome) - int(k):
        kmers[genome[i:i + int(k)]] = kmers.get(genome[i:i + int(k)], 0) + 1
        i += 1

    frequencies = []
    for kmer in sorted(kmers.keys()):
        frequencies.append(kmers[kmer])
    return frequencies


def pattern_to_number(pattern):
    """
    Assign a number to the pattern for its given lexigraphical position
    against all alphabetically sorted k-mers of the same length.
    :param pattern: k-mer pattern
    :return: lexigraphic position number
    """
    if not pattern:
        return 0
    return pattern_to_number(pattern[:-1]) * 4 + CONVERT[pattern[-1]]


def number_to_pattern(number, k):
    """
    Take a lexigraphical position number and return the k-mer it
    represents.
    :param number: lexigraphic position number
    :param k: length of k-mer
    :return: k-mer
    """
    if int(k) == 1:
        return PROFILE_INDEX[int(number)]
    return number_to_pattern(math.floor(int(number) / 4), int(k) - 1) + PROFILE_INDEX[int(number) % 4]


def get_kmer_counts(genome, k, d=0, rc=False):
    """
    Find counts for all approximate kmers
    :param genome: genome string
    :param k: k-mer size
    :param d: maximum distance value
    :param rc: toggle between whether reverse compliments count toward
        the final count
    :return: {approximate_kmer: count}
    """
    approximate_kmers = {}

    # Set initial kmer counts
    kmers = get_all_kmers(k)
    i = 0
    while i <= len(genome) - int(k):
        kmers[genome[i:i + int(k)]] = kmers.get(genome[i:i + int(k)], 0) + 1
        i += 1

    # Handle adding approximate counts if applicable
    for kmer in list(kmers.keys()):
        approximate_kmers[kmer] = approximate_kmers.get(kmer, 0) + kmers[kmer]  # Add count of kmer itself
        for other_kmer in [k for k in kmers.keys() if k != kmer]:
            # For efficiency only run hamming distances if either kmer has a count
            if kmers[kmer] or kmers[other_kmer]:
                # If hamming distance match, add the count of each to the others final count
                if hamming_distance(kmer, other_kmer) <= d:
                    approximate_kmers[kmer] = approximate_kmers.get(kmer, 0) + kmers[other_kmer]
                    approximate_kmers[other_kmer] = approximate_kmers.get(other_kmer, 0) + kmers[kmer]
        # For efficiency in not duplicating hamming distance calculations delete current key
        del kmers[kmer]

    # If reverse compliments is True add counts of the reverse compliment
    if rc:
        check = {}
        for kmer in approximate_kmers:
            # If the kmer is in check, it was already handled by its compliment
            if kmer in check:
                continue
            complement = reverse_complement(kmer)
            approximate_kmers[kmer] = approximate_kmers.get(kmer, 0) + approximate_kmers[complement]
            approximate_kmers[complement] = approximate_kmers[kmer]
            check[complement] = True

    return approximate_kmers


def get_frequent_kmers(genome, k, d=0, rc=False):
    """
    Find the most frequent k-mers for a genome genome.
    :param genome: genome string
    :param k: length of k-mer
    :param d: maximum hamming distance value if approximate method
    :param rc: reverse compliment processing
    :return: space delimited list of the most frequent kmers
    """
    kmers = get_kmer_counts(genome, int(k), int(d), rc)
    value = kmers[max(kmers, key=kmers.get)]
    return ' '.join(sorted([k for k, v in kmers.items() if v == value]))


def find_skew(genome):
    """
    Find the skew of a genome.
    :param genome: genome string
    :return: space delimited skew string
    """
    skew = [0]
    for base in genome:
        diff = 1 if base == 'G' else (-1 if base == 'C' else 0)
        skew.append(skew[-1] + diff)
    return ' '.join(str(i) for i in skew)


def find_min_skews(genome):
    """
    Find the minimum skew indices of a genome.
    :param genome: genome string
    :return: space delimited string of minimum skews
    """
    prev = 0
    cur_min = 0
    min_skews = []
    i = 0
    while i < len(genome):
        prev += 1 if genome[i] == 'G' else (-1 if genome[i] == 'C' else 0)
        if prev < cur_min:
            cur_min = prev
            min_skews = [i+1]
        elif prev == cur_min:
            min_skews.append(i+1)
        i += 1
    return ' '.join(str(i) for i in min_skews)


def find_max_skews(genome):
    """
    Find the maximum skew indices of a genome.
    :param genome: genome string
    :return: space delimited string of maximum skews
    """
    prev = 0
    cur_max = 0
    max_skews = []
    i = 0
    while i < len(genome):
        prev += 1 if genome[i] == 'G' else (-1 if genome[i] == 'C' else 0)
        if prev > cur_max:
            cur_max = prev
            max_skews = [i+1]
        elif prev == cur_max:
            max_skews.append(i+1)
        i += 1
    return ' '.join(str(i) for i in max_skews)


def motif_enumerations(genomes, k, d):
    """
    Find all kmers that occur within every genome in genomes with at
    most a hamming distance of d.
    :param genomes: list of genome strings
    :param k: length of k-mers
    :param d: maximum hamming distance
    :return: all kmers that occur in ever genome
    """
    kmers = {}  # records whether existing kmer matches are found in the current line
    patterns = {}  # lookup of observed patterns to kmers that can still satisfy being found in every line
    lookup = {}  # efficiency to recalculate strings at the end of each line by pointing which strings drop which kmers
    init = True  # toggle to differentiate initial action on first genome to actions on the following genomes

    for genome in genomes:

        # Explore the genome by pattern of length k for each index until the end
        i = 0
        while i < len(genome) - int(k) + 1:
            # If current pattern not in patterns find its matching kmers
            if genome[i:i+int(k)] not in patterns:
                # Find all potential mismatched kmers
                patterns[genome[i:i+int(k)]] = get_all_mismatched_kmers(genome[i:i+int(k)], int(d))

                # Setup lookup -> {kmer: {pattern: True}}
                if init:
                    # If this is the first line set all kmer options to the lookup
                    for kmer in patterns[genome[i:i+int(k)]]:
                        if kmer not in lookup:
                            lookup[kmer] = {}
                        lookup[kmer][genome[i:i+int(k)]] = True

                else:
                    # Only allow kmer options that already exist from previous lines
                    for kmer in list(patterns[genome[i:i+int(k)]].keys()):
                        if kmer in lookup:
                            lookup[kmer][genome[i:i+int(k)]] = True
                        # Remove kmers that have not appeared in every previous line from the pattern's lookup
                        else:
                            del patterns[genome[i:i+int(k)]][kmer]

            # For every kmer still matched to the pattern, set it to True for appearing in the line
            for kmer in patterns[genome[i:i+int(k)]]:
                kmers[kmer] = True

            i += 1

        # Handle end-of-line tasks
        for kmer in list(kmers.keys()):
            # If the kmer did not appear in this line remove it from all references
            if not kmers[kmer]:
                for pattern in lookup[kmer]:
                    del patterns[pattern][kmer]
                del lookup[kmer]
                del kmers[kmer]

            # Otherwise set it to False for processing on the next line
            else:
                kmers[kmer] = False

        init = False

    return ' '.join(sorted(kmers.keys()))


def median_string(genomes, k):
    """
    Find the k-mer with the minimum hamming distance score summed over
    all lines within the genomes.
    :param genomes: list of genome strings
    :param k: length of k-mers
    :return: motif k-mer with the lowest score the genome
    """
    kmers = {}  # records the score for all kmers
    patterns = {}  # lookup of observed patterns to kmers and their distance score
    lookup = {}  # store the minimum score for kmers in each line

    for genome in genomes:

        # Explore the genome by pattern of length k for each index until the end
        i = 0
        while i < len(genome) - int(k) + 1:
            # If current pattern not in patterns find its matching kmers
            if genome[i:i+int(k)] not in patterns:
                # Find all potential mismatched kmers
                patterns[genome[i:i+int(k)]] = get_all_mismatched_kmers(genome[i:i+int(k)], int(k))

            # If this is the first pattern for the genome line the lookup is equal to the pattern scores
            if not i:
                lookup = copy.deepcopy(patterns[genome[i:i+int(k)]])

            else:
                # Check to see if a new minimum applies to each kmer
                for kmer in patterns[genome[i:i+int(k)]]:
                    if patterns[genome[i:i+int(k)]][kmer] < lookup[kmer]:
                        lookup[kmer] = patterns[genome[i:i+int(k)]][kmer]

            i += 1

        # Handle end-of-line summing
        for kmer in lookup:
            kmers[kmer] = kmers.get(kmer, 0) + lookup[kmer]

    return min(kmers, key=kmers.get)


def get_profile_dict(profile_matrix):
    """
    Take a list of space delimited probabilities and convert it to a
    dictionary of {column: {base: probability}}.
    :param profile_matrix: list of space delimited probabilities
    :return: dictionary representation of profile matrix
    """
    profile = {}

    # Convert space delimited entries to actual lists
    matrix = []
    for row in profile_matrix:
        matrix.append([float(f) for f in row.split()])

    # Process by column
    j = 0
    while j < len(matrix[0]):

        # Add column to profile
        profile[j] = {}

        # Add each base and its probability
        i = 0
        while i < len(matrix):
            profile[j][PROFILE_INDEX[i]] = matrix[i][j]
            i += 1

        j += 1

    return profile


def greedy_motif_search(genomes, k, laplace=True):
    """
    Greedy algorithm solution to motif finding; find T motifs that
    produce the lowest score for the genome.
    :param genomes: set of genome strings
    :param k: length of kmers
    :param laplace: smoothing parameter for profiles
    :return: set of motifs with the lowest score
    """
    # Build initial best motifs from first kmers of each line
    best_motifs = dict((i, {}) for i in range(0, int(k)))
    i = 0
    while i < len(genomes):
        j = 0
        while j < int(k):
            best_motifs[j][i] = genomes[i][j]
            j += 1
        i += 1
    best_score = score_matrix(best_motifs)

    # Iterate through each kmer in the first genome string
    j = 0
    while j < len(genomes[0]) - int(k) + 1:
        # Set motif matrix to the first genome string's kmer
        motifs = dict((i, {0: genomes[0][j+i]}) for i in range(0, int(k)))

        # Iterate through each following genome string and update the motif and profile matrices
        i = 1
        while i < len(genomes):
            new_motif = most_probable_in_profile(genomes[i], int(k), build_profile(motifs, laplace))
            for x in range(0, int(k)):
                motifs[x][i] = new_motif[x]
            i += 1

        # If the score for the motifs is better than the best score set motifs as the best motifs
        if score_matrix(motifs) < best_score:
            best_motifs = motifs
            best_score = score_matrix(motifs)

        j += 1

    return convert_matrix_to_tuple(best_motifs)


def randomized_motif_search(genomes, k, laplace=True, trials=1000):
    """
    Randomized algorithm solution to motif finding; find T motifs that
    produce the lowest score for the genome.
    :param genomes: set of genome strings
    :param k: length of kmers
    :param laplace: smoothing parameter for profiles
    :param trials: number of times the random search will be completed
    :return: set of motifs with the lowest score
    """
    best_score = sys.maxsize
    best_motifs = {}

    # Process each trial
    for t in range(0, int(trials)):

        # Build initial best motifs from random kmers of each line
        motifs = dict((i, {}) for i in range(0, int(k)))
        i = 0
        while i < len(genomes):
            r = random.randint(0, len(genomes[i]) - int(k))
            j = 0
            while j < int(k):
                motifs[j][i] = genomes[i][r+j]
                j += 1
            i += 1

        # Iteratively find the most probable motifs given the previous profile until the score cannot be improved
        while motifs:
            profile = build_profile(motifs, laplace)

            # Iterate through each following genome string and update the motif
            i = 0
            new_motifs = dict((i, {}) for i in range(0, int(k)))
            while i < len(genomes):
                new_motif = most_probable_in_profile(genomes[i], int(k), profile)
                for x in range(0, int(k)):
                    new_motifs[x][i] = new_motif[x]
                i += 1

            # If the score for the motifs is better than the best score set motifs as the best motifs
            if score_matrix(new_motifs) < score_matrix(motifs):
                motifs = new_motifs

            # A better motif matrix was not found so terminate this trial
            else:
                # If this trial did better than all other previous trials record motifs and score
                if score_matrix(motifs) < best_score:
                    best_motifs = motifs
                    best_score = score_matrix(motifs)
                motifs = None

    return convert_matrix_to_tuple(best_motifs)


def randomized_gibbs_motif_search(genomes, k, gibbs=100, laplace=True, trials=1000):
    """
    Randomized algorithm solution to motif finding; find T motifs that
    produce the lowest score for the genome.
    :param genomes: set of genome strings
    :param k: length of kmers
    :param gibbs: number of times the gibbs search is run
    :param laplace: smoothing parameter for profiles
    :param trials: number of times the random search will be completed
    :return: set of motifs with the lowest score
    """
    best_score = sys.maxsize
    best_motifs = {}

    # Process each trial
    for t in range(0, int(trials)):

        # Build initial best motifs from random kmers of each line
        motifs = dict((i, {}) for i in range(0, int(k)))
        i = 0
        while i < len(genomes):
            r = random.randint(0, len(genomes[i]) - int(k))
            j = 0
            while j < int(k):
                motifs[j][i] = genomes[i][r+j]
                j += 1
            i += 1

        # Replace one motif at a time by random selection for gibbs number of times
        for g in range(0, int(gibbs)):
            # Select a random line for which to replace using an updated profile matrix
            r = random.randint(0, len(genomes) - 1)

            # Delete existing kmer motif of randomly selected line
            for col in motifs:
                del motifs[col][r]

            # Build a profile from the motifs and select a new best motif from the randomly selected line
            new_motif = random_probable_in_profile(genomes[r], int(k), build_profile(motifs, laplace))
            for x in range(0, int(k)):
                motifs[x][r] = new_motif[x]

            # If the score for the motifs is better than the best score set motifs as the best motifs
            if score_matrix(motifs) < best_score:
                best_motifs = copy.deepcopy(motifs)
                best_score = score_matrix(motifs)

    return convert_matrix_to_tuple(best_motifs)


def most_probable_in_profile(genome, k, profile):
    """
    Find the most probably motif given a profile matrix and a length of
    k within the genome.
    :param genome: genome string
    :param k: length of kmers/motif
    :param profile_matrix: profile probabilities for each base
    :return: most probable motif
    """
    max_prob = (0, genome[:int(k)])

    # Explore the genome by kmer of length k for each index until the end
    i = 0
    while i < len(genome) - int(k) + 1:

        # Consider the kmer by base, stop if the current sum is less than the max prob
        j = 0
        prob = 1
        while j < int(k):

            # Multiply the prob by the profile probability of the current j index and observed base
            try:
                prob *= profile[j][genome[i+j]]
            except KeyError:
                prob = 0

            # If prob probability is less than the max_prob, stop processing kmer
            if prob < max_prob[0]:
                j = k
            j += 1

        # If prob exceeded the max_prob replace with the current kmer and its probability
        if prob > max_prob[0]:
            max_prob = (prob, genome[i:i+int(k)])
        i += 1

    return max_prob[1]


def random_probable_in_profile(genome, k, profile):
    """
    Find a random motif given probabilities of a profile matrix and a
    length of k within the genome.
    :param genome: genome string
    :param k: length of kmers/motif
    :param profile_matrix: profile probabilities for each base
    :return: random probable motif
    """
    kmers = []
    probs = []
    total = 0.0

    # Explore the genome by kmer of length k for each index until the end
    i = 0
    while i < len(genome) - int(k) + 1:

        # Consider the kmer by base, stop if the current sum is less than the max prob
        j = 0
        prob = 1
        while j < int(k):

            # Multiply the prob by the profile probability of the current j index and observed base
            try:
                prob *= profile[j][genome[i+j]]
            except KeyError:
                prob = 0

            j += 1

        # Add kmer and probability to lists
        kmers.append(genome[i:i+k])
        probs.append(prob)
        total += prob
        i += 1

    # Normalize probabilities
    normalized_probs = []
    for prob in probs:
        normalized_probs.append(prob/total)

    return choice(kmers, p=normalized_probs)


def build_profile(motif_matrix, laplace=True):
    """
    Convert a motif_matrix to a profile matrix dict with probabilities.
    :param motif_matrix: a motif matrix in dictionary form
    :param laplace: smoothing parameter for profiles
    :return: profile matrix in dictionary form
    """
    profile = {}
    for col in motif_matrix:

        # Collect the distribution of bases for the column
        distr = dict((b, 1) for b in BASES) if laplace else {}
        for row in motif_matrix[col]:
            distr[motif_matrix[col][row]] = distr.get(motif_matrix[col][row], 0) + 1

        # Convert the distributions into probabilities
        for base in distr:
            distr[base] = distr.get(base, 0) / (len(motif_matrix[col]) + 4 if laplace else len(motif_matrix[col]))

        profile[col] = copy.deepcopy(distr)

    return profile


def score_matrix(motif_matrix):
    """
    Take a motif matrix in dictionary form {column: {row: base}} and
    output the score for the motif matrix.
    :param motif_matrix: a motif matrix in dictionary form
    :return: score for the motif matrix
    """
    score = 0
    for col in motif_matrix:

        # Collect the distribution of bases for the column
        distr = {}
        for row in motif_matrix[col]:
            distr[motif_matrix[col][row]] = distr.get(motif_matrix[col][row], 0) + 1

        # Score by taking the total rows less the count of the most frequent base
        score += (len(motif_matrix[col]) - distr[max(distr, key=distr.get)])
    return score


def convert_matrix_to_tuple(motif_matrix):
    """
    Take a motif matrix in dictionary form {column: {row: base}} and
    output a tuple of the motifs
    :param motif_matrix: a motif matrix in dictionary form
    :return: tuple of motifs
    """
    motifs = {}
    for col in sorted(motif_matrix.keys()):
        for row in motif_matrix[col]:
            motifs[row] = motifs.get(row, '') + motif_matrix[col][row]
    return tuple(v for k, v in sorted(motifs.items()))


def distance_between_pattern_and_strings(pattern, genomes):
    """
    Find the total distance between a pattern and each string in
    genomes. Take the minimum distance for each string. This is NOT
    used in median string because it is very inefficient.
    :param pattern: k-mer pattern
    :param genomes: set of genomes
    :return: distance
    """
    score = 0
    for genome in genomes:
        i = 0
        current = sys.maxsize
        while i <= len(genome) - len(pattern):
            if hamming_distance(pattern, genome[i:i+len(pattern)]) < current:
                current = hamming_distance(pattern, genome[i:i+len(pattern)])
            i += 1
        score += current
    return score


lines = sys.stdin.read().splitlines()
print(number_to_pattern(*lines))
