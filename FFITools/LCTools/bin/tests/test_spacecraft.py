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
Test tsig.spacecraft API in a general way.
"""

import os
import sys

from tsig.spacecraft import Spacecraft

from unittest import TestCase

try:
    from test import test_support
except ImportError:
    from test import support as test_support

class SpacecraftApiTestCase(TestCase):
    """Unit tests for all classes in the tsig.spacecraft module."""

    def test_01_default_config(self):
        """read spacecraft configuration"""
        sc = Spacecraft()
        self.assertEqual(sc.camera[0].label, "camera_1")
        self.assertEqual(sc.camera[1].label, "camera_2")
        self.assertEqual(sc.camera[2].label, "camera_3")
        self.assertEqual(sc.camera[3].label, "camera_4")
        self.assertEqual(sc.camera[0].serial_number, "sn05")
        self.assertEqual(sc.camera[1].serial_number, "sn06")
        self.assertEqual(sc.camera[2].serial_number, "sn07")
        self.assertEqual(sc.camera[3].serial_number, "sn09")
        self.assertEqual(sc.camera[0].ccd[0].label, "ccd_1")
        self.assertEqual(sc.camera[0].offset_angle, -36)
        self.assertEqual(sc.camera[1].offset_angle, -12)
        self.assertEqual(sc.camera[2].offset_angle, 12)
        self.assertEqual(sc.camera[3].offset_angle, 36)
        self.assertEqual(sc.camera[0].ccd[0].rows, 2058)
        self.assertEqual(sc.camera[0].ccd[0].cols, 2048)

    def test_02_spacecraft_pointing(self):
        """point a spacecraft, check the cameras"""
        spacecraft = Spacecraft()
        spacecraft.set_pointing(0, 0, 0)
        ra, dec, roll = spacecraft.camera[0].get_pointing()
        self.assertEqual([ra, dec, roll], [0.0, -36.0, 0.0])
        ra, dec, roll = spacecraft.camera[1].get_pointing()
        self.assertEqual([ra, dec, roll], [0.0, -12.0, 0.0])
        ra, dec, roll = spacecraft.camera[2].get_pointing()
        self.assertEqual([ra, dec, roll], [0.0, 12.0, 0.0])
        ra, dec, roll = spacecraft.camera[3].get_pointing()
        self.assertEqual([ra, dec, roll], [0.0, 36.0, 0.0])


if __name__ == '__main__':
    test_support.run_unittest(SpacecraftApiTestCase)
