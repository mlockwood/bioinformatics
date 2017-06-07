#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


class KMeansCluster(object):

    def __init__(self, data, k, set_start):
        self.data = data
        self.k = int(k)
        self.distances = {}
        self.centers = {*[set_start]}  # for proper random use {*[next(iter(data))]}
        self.farthest_first_traversal()

    @staticmethod
    def load_data_from_strings(lines):
        """
        Intake lines of strings with space separated values and upload them
        to the dataset.
        :param lines: string delimited values
        :return: set of data points
        """
        data = set()
        for line in lines:
            data.add(tuple(float(f) for f in line.rstrip().split()))
        return data

    @staticmethod
    def euclidean_distance(v, w):
        """
        Take two vectors and return their squared euclidean distance. The
        square root is not taken to save computation time and does not
        impact results for the k-means clustering.
        :param v: first vector (list or tuple)
        :param w: second vector (list or tuple)
        :return: euclidean distance
        """
        distance = 0
        i = 0
        while i < len(v):
            distance += (v[i] - w[i]) ** 2
            i += 1
        return distance

    def farthest_first_traversal(self):
        """
        Initialize the centers for the k means clustering function.
        :return: set of centers
        """
        last_center = next(iter(self.centers))
        while len(self.centers) < self.k:

            # Add the farthest point from any center as a new center
            farthest = (0, None)
            for point in self.data:

                # Store lookup of euclidean distance between point and centers
                if point not in self.distances:
                    self.distances[point] = {}
                self.distances[point][last_center] = KMeansCluster.euclidean_distance(point, last_center)

                # Compare minimum center distance calculation to current lead
                closest = self.distances[point][min(self.distances[point], key=self.distances[point].get)]
                if closest > farthest[0]:
                    farthest = (closest, point)

            self.centers.add(farthest[1])
            last_center = farthest[1]

    def print_ready(self):
        out = ''
        for center in self.centers:
            out += ' '.join(str(s) for s in center)
            out += '\n'
        return out


lines = sys.stdin.read().splitlines()
data = KMeansCluster.load_data_from_strings(lines[1:])
print(KMeansCluster(data, int(lines[0].split()[0]), tuple(float(f) for f in lines[1].rstrip().split())).print_ready())
