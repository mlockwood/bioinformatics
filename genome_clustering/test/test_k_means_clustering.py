#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from genome_clustering.k_means_clustering import *


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


data_lines = ['0.0 0.0', '5.0 5.0', '0.0 5.0', '1.0 1.0', '2.0 2.0', '3.0 3.0', '1.0 2.0']
kmeans_obj = KMeansCluster(KMeansCluster.load_data_from_strings(data_lines), 3)


def test_euclidean_distance():
    assert KMeansCluster.euclidean_distance((1, 2), (3, 2)) == 4


def test_farthest_first_traversal():
    assert (5, 5) in KMeansCluster(KMeansCluster.load_data_from_strings(data_lines), 3).centers
