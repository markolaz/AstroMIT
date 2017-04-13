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
Test tsig.camera API in a general way.
"""

import os
import sys

from unittest import TestCase

try:
    from test import test_support
except ImportError:
    from test import support as test_support

class CameraApiTestCase(TestCase):
    """Unit tests for all classes in the tsig.camera module.
    
    We test that the importing is working and that the basic setup init
    functions are all working, but no more. Further testing for each class
    should be done in it's own test file.
    """
    def test_01_camera_config(self):
        from tsig.camera import Camera
        from tsig.camera.camera import Camera
        camera = Camera('test')
        self.assertEqual(camera.label, "test")
        self.assertEqual(camera.serial_number, 0)

    def test_02_camera_pointing(self):
        from tsig.camera import Camera
        from tsig.camera.camera import Camera
        camera = Camera('test')
        camera.point_camera(10, 20, 30)
        pt = camera.get_pointing()
        self.assertEqual(str(pt), "(10, 20, 30)")

    def test_03_ccd_config(self):
        from tsig.camera import CCD
        from tsig.camera.ccd import CCD
        ccd = CCD('test')
        self.assertEqual(ccd.label, "test")
        self.assertEqual(ccd.rows, 2058)
        self.assertEqual(ccd.cols, 2048)


if __name__ == '__main__':
    test_support.run_unittest(CameraApiTestCase)
