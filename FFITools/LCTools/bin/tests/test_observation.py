#!/usr/bin/env python
#
# Copyright (C) 2017 - Massachusetts Institute of Technology (MIT)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Tests for the tsig Observation object.
"""

import os
import sys

from unittest import TestCase

try:
    from test import test_support
except ImportError:
    from test import support as test_support

class ObservationTestCase(TestCase):
    def test_01_import(self):
        from tsig.observation import Observation
        obs = Observation()


if __name__ == '__main__':
    test_support.run_unittest(ObservationTestCase)
