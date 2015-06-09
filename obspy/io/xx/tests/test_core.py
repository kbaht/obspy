#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA @UnusedWildImport

import os
import unittest

import numpy as np

from obspy import UTCDateTime, read
from obspy.io.xx.core import _is_xx, _read_xx


class CoreTestCase(unittest.TestCase):
    """
    Baykal XX file test suite.
    """

    def setUp(self):
        # Directory where the test files are located
        self.path = os.path.join(os.path.dirname(__file__), 'data')

    def test_is_xx(self):
        """
        Testing XX file format.
        """
