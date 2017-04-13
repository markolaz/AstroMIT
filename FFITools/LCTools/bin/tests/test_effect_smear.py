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
Test tsig.effects.smear
"""

import os
import sys
import numpy as np
from tsig.effects.smear import Smear

from base import TsigTestCase, test_support

class SmearTestCase(TsigTestCase):
    def test_01_smear(self):
        smear = Smear()

    def test_02_smear_square(self):
        """ensure that smear works on square images"""
        image = np.zeros(shape=(100, 100), dtype=np.int32)
        smear = Smear()
        result = smear.apply_to_image(image)
        self.assertEqual(str(result), str(image))

    def test_02_smear_non_square(self):
        """ensure that smear works on non-square images"""
        image = np.zeros(shape=(100, 200), dtype=np.int32)
        smear = Smear()
        result = smear.apply_to_image(image)
        self.assertEqual(str(result), str(image))


if __name__ == '__main__':
    test_support.run_unittest(SmearTestCase)

