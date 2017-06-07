#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from genome_clustering.k_means_clustering import *


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def test_euclidean_distance():
    assert KMeansCluster.euclidean_distance((1, 2), (3, 2)) == 4


def test_squared_error_distortion():
    center_lines = ['2.31 4.55', '5.96 9.08']
    data_lines = ['3.42 6.03', '6.23 8.25', '4.76 1.64', '4.47 4.33', '3.95 7.61', '8.93 2.97', '9.74 4.03',
                  '1.73 1.28', '9.72 5.01', '7.27 3.77']
    k_means_obj = KMeansCluster(KMeansCluster.load_points_from_strings(data_lines),
                                2,
                                KMeansCluster.load_points_from_strings(center_lines))
    assert 18.245 <= k_means_obj.squared_error_distortion() <= 18.246


def test_lloyd_algorithm():
    data_lines = ['1.3 1.1', '1.3 0.2', '0.6 2.8', '3.0 3.2', '1.2 0.7', '1.4 1.6', '1.2 1.0', '1.2 1.1', '0.6 1.5',
                  '1.8 2.6', '1.2 1.3', '1.2 1.0', '0.0 1.9']
    k_means_obj = KMeansCluster(KMeansCluster.load_points_from_strings(data_lines), 2)
    assert (1.800, 2.867) in k_means_obj.centers and (1.060, 1.140) in k_means_obj.centers
