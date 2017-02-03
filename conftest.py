#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true", help="run slow tests")
    parser.addoption("--random", action="store_true", help="run tests with random or inconsistent outputs")