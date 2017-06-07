#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import copy
import sys


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


class KMeansCluster(object):

    def __init__(self, data, k, centers=None):
        self.data = data
        self.k = int(k)
        self.distances = {}
        self.centers = {*[next(iter(data))]} if not centers else centers
        if not centers:
            self.lloyd_algorithm()

    @staticmethod
    def load_points_from_strings(lines):
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

    @staticmethod
    def center_of_gravity(points):
        """
        Find the center of gravity for a series of points.
        :param points: points that form a cluster
        :return: the point/center of gravity
        """
        print(points)
        # First sort point information by dimension
        gravity = dict((i, []) for i in range(len(points[0])))
        for point in points:
            i = 0
            while i < len(point):
                gravity[i] = gravity.get(i) + [point[i]]
                i += 1

        # Average each dimension and then return the point
        point = []
        for dimension in sorted(gravity.keys()):
            point.append(sum(gravity[dimension]) / len(gravity[dimension]))

        print(point)
        return tuple(point)

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

    def squared_error_distortion(self):
        """
        Find the squared error distortion of the current points and
        centers.
        :return: The squared error distortion 
        """
        distortion = 0
        for point in self.data:
            # If centers were provided calculate the distances
            # Store lookup of euclidean distance between point and centers
            for center in self.centers:
                if point not in self.distances:
                    self.distances[point] = {}
                self.distances[point][center] = KMeansCluster.euclidean_distance(point, center)

            # Use minimum center distance to find squared distortion
            distortion += self.distances[point][min(self.distances[point], key=self.distances[point].get)]
        return distortion / len(self.data)  # remember to return an average

    def lloyd_algorithm(self):
        self.farthest_first_traversal()
        prev_centers = None
        while prev_centers != self.centers:
            prev_centers = copy.deepcopy(self.centers)

            # Assign points to a center
            grouped_centers = dict((center, []) for center in self.centers)
            for point in self.data:
                self.distances[point] = {}
                # Store lookup of euclidean distance between point and centers
                for center in self.centers:
                    self.distances[point][center] = KMeansCluster.euclidean_distance(point, center)

                # Find closest center
                grouped_centers[min(self.distances[point], key=self.distances[point].get)] = grouped_centers.get(
                    min(self.distances[point], key=self.distances[point].get)) + [point]

            # Find new centers
            self.centers = set()
            for center in grouped_centers:
                self.centers.add(KMeansCluster.center_of_gravity(grouped_centers[center]))

    def print_centers(self):
        out = ''
        for center in self.centers:
            out += ' '.join(center)
            out += '\n'
        return out


# lines = sys.stdin.read().splitlines()
# k = int(lines[0].split()[0])
# data = KMeansCluster.load_points_from_strings(lines[1:])
# print(KMeansCluster(data, k).print_centers())
