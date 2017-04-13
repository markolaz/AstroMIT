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
Test transformations.
"""

import os
import sys
import numpy as np

from unittest import TestCase

try:
    from test import test_support
except ImportError:
    from test import support as test_support

from tsig.spacecraft.transformer import Transformer

class TransformationsApiTestCase(TestCase):
    """Unit tests for all classes in the tsig.transformer."""

    def get_cam_geometry(self, camera_number):
        """sample camera geometry for testing purposes"""
        return {
            'angle': {1: -36.0, 2: -12.0, 3: 12.0, 4: 36.0}.get(camera_number),
            'number': camera_number,
            'fov_radius': 16.970562748477143,
            'focal_length': 146.985,
            'a_0':   1.00000140,
            'a_2':   0.28174612,
            'a_4':  -0.59667259,
            'a_6':   9.17151267,
            'a_8':  -4.36928235,
            'b_0':   0.99999822,
            'b_2':  -0.28159976,
            'b_4':   0.82104074,
            'b_6': -10.22208978,
            'b_8':  22.83962334}

    def get_ccd_geometries(self, ccd_number):
        """sample ccd geometries for testing purposes"""
        return {
            1: {'rows': 2048,
                'cols': 2048,
                'rotation': 180.0,
                'x_0': 30.82, # mm
                'y_0': 30.82, # mm
                'pixel_size_x': 0.015, # mm
                'pixel_size_y': 0.015}, # mm
            2: {
                'rows': 2048,
                'cols': 2048,
                'rotation': 180.0,
                'x_0': 0.0,
                'y_0': 30.82,
                'pixel_size_x': 0.015,
                'pixel_size_y': 0.015},
            3: {
                'rows': 2048,
                'cols': 2048,
                'rotation': 0.0,
                'x_0': -30.82,
                'y_0': -30.82,
                'pixel_size_x': 0.015,
                'pixel_size_y': 0.015},
            4: {
                'rows': 2048,
                'cols': 2048,
                'rotation': 0.0,
                'x_0': 0.0,
                'y_0': -30.82,
                'pixel_size_x': 0.015,
                'pixel_size_y': 0.015}
        }

    def test_01_eq_to_sc_matrix(self):
        xf = Transformer()
        m = xf.get_eq_sc_matrix(0, 0, 0)
        self.assertEqual(str(m),
                         "[[ 0. -0.  1.]\n"
                         " [ 0. -1. -0.]\n"
                         " [ 1.  0.  0.]]")
        m = xf.get_eq_sc_matrix(45, 0, 0)
        self.assertEqual(str(m),
                         "[[ 0.         -0.          1.        ]\n"
                         " [ 0.70710678 -0.70710678 -0.        ]\n"
                         " [ 0.70710678  0.70710678  0.        ]]")
        m = xf.get_eq_sc_matrix(0, 45, 0)
        self.assertEqual(str(m),
                         "[[-0.70710678 -0.          0.70710678]\n"
                         " [ 0.         -1.         -0.        ]\n"
                         " [ 0.70710678  0.          0.70710678]]")
        m = xf.get_eq_sc_matrix(0, 0, 45)
        self.assertEqual(str(m),
                         "[[ 0.         -0.70710678  0.70710678]\n"
                         " [ 0.         -0.70710678 -0.70710678]\n"
                         " [ 1.          0.          0.        ]]")

    def test_02_sc_to_cam_matrix(self):
        xf = Transformer()
        m = xf.get_sc_cam_matrix(-36.0)
        self.assertEqual(str(m),
                         "[[ 0.         -1.          0.        ]\n"
                         " [ 0.80901699  0.          0.58778525]\n"
                         " [-0.58778525  0.          0.80901699]]")
        m = xf.get_sc_cam_matrix(-12.0)
        self.assertEqual(str(m),
                         "[[ 0.         -1.          0.        ]\n"
                         " [ 0.9781476   0.          0.20791169]\n"
                         " [-0.20791169  0.          0.9781476 ]]")
        m = xf.get_sc_cam_matrix(12.0)
        self.assertEqual(str(m),
                         "[[ 0.         -1.          0.        ]\n"
                         " [ 0.9781476   0.         -0.20791169]\n"
                         " [ 0.20791169  0.          0.9781476 ]]")
        m = xf.get_sc_cam_matrix(36.0)
        self.assertEqual(str(m),
                         "[[ 0.         -1.          0.        ]\n"
                         " [ 0.80901699  0.         -0.58778525]\n"
                         " [ 0.58778525  0.          0.80901699]]")

    def test_10_spacecraft_to_camera_pointing(self):
        # point the spacecraft then check each camera pointing
        xf = Transformer()
        m = xf.spacecraft_to_camera_pointing(0, 0, 0, -36.0)
        self.assertEqual(m, (0.0, -36.0, 0.0))
        m = xf.spacecraft_to_camera_pointing(0, 0, 0, -12.0)
        self.assertEqual(m, (0.0, -12.0, 0.0))
        m = xf.spacecraft_to_camera_pointing(0, 0, 0, 12.0)
        self.assertEqual(m, (0.0, 12.0, 0.0))
        m = xf.spacecraft_to_camera_pointing(0, 0, 0, 36.0)
        self.assertEqual(m, (0.0, 36.0, 0.0))

        m = xf.spacecraft_to_camera_pointing(
            147.418172456,
            -46.0280296711,
            -214.75909812,
            -36.0)
        print m

    def test_11_celestial_to_camera(self):
        xf = Transformer()
        # FIXME: verify the field angles
        ra, dec, fax, fay = xf.celestial_to_camera(0, 0, 0, 0, 0, -36.0)
        self.assertEqual(str([ra, dec]), "[array([ 90.]), array([ 54.])]")
        self.assertEqual(str([fax, fay]), "[array([ 0.]), array([ 36.])]")
        ra, dec, fax, fay = xf.celestial_to_camera(0, 0, 0, 0, 0, -12.0)
        self.assertEqual(str([ra, dec]), "[array([ 90.]), array([ 78.])]")
        self.assertEqual(str([fax, fay]), "[array([ 0.]), array([ 12.])]")
        ra, dec, fax, fay = xf.celestial_to_camera(0, 0, 0, 0, 0, 12.0)
        self.assertEqual(str([ra, dec]), "[array([-90.]), array([ 78.])]")
        self.assertEqual(str([fax, fay]), "[array([ 0.]), array([-12.])]")
        ra, dec, fax, fay = xf.celestial_to_camera(0, 0, 0, 0, 0, 36.0)
        self.assertEqual(str([ra, dec]), "[array([-90.]), array([ 54.])]")
        self.assertEqual(str([fax, fay]), "[array([ 0.]), array([-36.])]")

        # a star near orion at 72.46, 6.961
        ra, dec, fax, fay = xf.celestial_to_camera(
            72.46, 6.961, 83.0, 12.0, 0, -36.0)
        self.assertEqual(str([ra, dec]),
                         "[array([ 109.68111887]), array([ 57.37502272])]")
        ra, dec, fax, fay = xf.celestial_to_camera(
            72.46, 6.961, 83.0, 12.0, 0, -12.0)
        self.assertEqual(str([ra, dec]),
                         "[array([ 146.27834281]), array([ 77.39052882])]")
        ra, dec, fax, fay = xf.celestial_to_camera(
            72.46, 6.961, 83.0, 12.0, 0, 12.0)
        self.assertEqual(str([ra, dec]),
                         "[array([-122.39127858]), array([ 70.18723765])]")
        ra, dec, fax, fay = xf.celestial_to_camera(
            72.46, 6.961, 83.0, 12.0, 0, 36.0)
        self.assertEqual(str([ra, dec]),
                         "[array([-105.74268186]), array([ 47.9924046])]")

    def test_11_celestial_to_camera_vector(self):
        xf = Transformer()
        ra, dec, fax, fay = xf.celestial_to_camera(
            [0, 10], [0, 20], 0, 0, 0, -36.0)
        self.assertEqual(str([ra, dec]),
                         "[array([ 90.       ,  78.7540829]),"
                         " array([ 54.        ,  33.20548617])]")
        ra, dec, fax, fay = xf.celestial_to_camera(
            [0, 10], [0, 20], 45, 90, 30, -36.0)
        self.assertEqual(str([ra, dec]),
                         "[array([-167.76892409, -171.96306913]),"
                         " array([  8.75038443,  30.67235662])]")

    def test_12_camera_to_focal(self):
        xf = Transformer()
        x, y = xf.camera_to_focal(0, 78, self.get_cam_geometry(1))
        self.assertEqual(str([x, y]), "[-31.628174060007204, -0.0]")
        x, y = xf.camera_to_focal(45, 78, self.get_cam_geometry(1))
        self.assertEqual(str([x, y]),
                         "[-22.364496354379554, -22.364496354379551]")
        x, y = xf.camera_to_focal(90, 78, self.get_cam_geometry(1))
        self.assertEqual(str([x, y]),
                         "[-1.9366671062731585e-15, -31.628174060007204]")
        x, y = xf.camera_to_focal(90, 54, self.get_cam_geometry(1))
        self.assertEqual(str([x, y]),
                         "[-1.3027233534502548e-14, -212.75086896193443]")

    def test_12_camera_to_focal_vector(self):
        xf = Transformer()
        x, y = xf.camera_to_focal([0, 45, 90, 90], [78, 78, 78, 54],
                                  self.get_cam_geometry(1))
        self.assertEqual(
            str([x, y]),
            "[array([ -3.16281741e+01,  -2.23644964e+01,  -1.93666711e-15,\n"
            "        -1.30272335e-14]),"
            " array([  -0.        ,  -22.36449635,  -31.62817406,"
            " -212.75086896])]")

    def test_13_celestial_to_focal(self):
        xf = Transformer()
        x, y, fax, fay = xf.celestial_to_focal(
            0, 0, 0, 0, 0, self.get_cam_geometry(1))
        self.assertEqual(str([x, y]),
                         "[array([ -1.30272335e-14]), array([-212.75086896])]")
        x, y, fax, fay = xf.celestial_to_focal(
            0, 0, 0, 0, 0, self.get_cam_geometry(2))
        self.assertEqual(str([x, y]),
                         "[array([ -1.93666711e-15]), array([-31.62817406])]")
        x, y, fax, fay = xf.celestial_to_focal(
            0, 0, 0, 0, 0, self.get_cam_geometry(3))
        self.assertEqual(str([x, y]),
                         "[array([ -1.93666711e-15]), array([ 31.62817406])]")
        x, y, fax, fay = xf.celestial_to_focal(
            0, 0, 0, 0, 0, self.get_cam_geometry(4))
        self.assertEqual(str([x, y]),
                         "[array([ -1.30272335e-14]), array([ 212.75086896])]")

    def test_13_celestial_to_focal_vector(self):
        xf = Transformer()
        x, y, fax, fay = xf.celestial_to_focal([0, 10, 20], [0, 20, 10], 0, 0, 0,
                                     self.get_cam_geometry(4))
        self.assertEqual(
            str([x, y]),
            "[array([ -1.30272335e-14,  -2.61319658e+01,  -8.36703519e+01]),"
            " array([ 212.75086896,   42.79838937,  100.22382439])]")

    def test_14_focal_to_pixel_cam_1(self):
        xf = Transformer()
        x, y, ccd_id = xf.focal_to_pixel(0.001, 0.001, 1,
                                              self.get_ccd_geometries(1))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([2056]), array([2056]), array([3])]")
        x, y, ccd_id = xf.focal_to_pixel(0.001, -0.001, 1,
                                              self.get_ccd_geometries(1))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([1]), array([2056]), array([2])]")
        x, y, ccd_id = xf.focal_to_pixel(-0.001, 0.001, 1,
                                              self.get_ccd_geometries(1))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([1]), array([2056]), array([4])]")
        x, y, ccd_id = xf.focal_to_pixel(-0.001, -0.001, 1,
                                              self.get_ccd_geometries(1))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([2056]), array([2056]), array([1])]")

    def test_14_focal_to_pixel_cam_2(self):
        xf = Transformer()
        x, y, ccd_id = xf.focal_to_pixel(0.001, 0.001, 2,
                                              self.get_ccd_geometries(2))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([2056]), array([2056]), array([3])]")
        x, y, ccd_id = xf.focal_to_pixel(0.001, -0.001, 2,
                                              self.get_ccd_geometries(2))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([1]), array([2056]), array([2])]")
        x, y, ccd_id = xf.focal_to_pixel(-0.001, 0.001, 2,
                                              self.get_ccd_geometries(2))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([1]), array([2056]), array([4])]")
        x, y, ccd_id = xf.focal_to_pixel(-0.001, -0.001, 2,
                                              self.get_ccd_geometries(2))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([2056]), array([2056]), array([1])]")

    def test_14_focal_to_pixel_cam_3(self):
        xf = Transformer()
        x, y, ccd_id = xf.focal_to_pixel(0.001, 0.001, 3,
                                              self.get_ccd_geometries(3))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([2056]), array([2056]), array([1])]")
        x, y, ccd_id = xf.focal_to_pixel(0.001, -0.001, 3,
                                              self.get_ccd_geometries(3))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([1]), array([2056]), array([4])]")
        x, y, ccd_id = xf.focal_to_pixel(-0.001, 0.001, 3,
                                              self.get_ccd_geometries(3))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([1]), array([2056]), array([2])]")
        x, y, ccd_id = xf.focal_to_pixel(-0.001, -0.001, 3,
                                              self.get_ccd_geometries(3))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([2056]), array([2056]), array([3])]")

    def test_14_focal_to_pixel_cam_4(self):
        xf = Transformer()
        x, y, ccd_id = xf.focal_to_pixel(0.001, 0.001, 4,
                                              self.get_ccd_geometries(4))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([2056]), array([2056]), array([1])]")
        x, y, ccd_id = xf.focal_to_pixel(0.001, -0.001, 4,
                                              self.get_ccd_geometries(4))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([1]), array([2056]), array([4])]")
        x, y, ccd_id = xf.focal_to_pixel(-0.001, 0.001, 4,
                                              self.get_ccd_geometries(4))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([1]), array([2056]), array([2])]")
        x, y, ccd_id = xf.focal_to_pixel(-0.001, -0.001, 4,
                                              self.get_ccd_geometries(4))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([2056]), array([2056]), array([3])]")

#    def test_14_focal_to_pixel_vector(self):
#        x, y, ccd_id = xf.focal_to_pixel([1, 2, 3], [1, 2, 3], 1)
#        self.assertEqual(str([x, y, ccd_id]), "")
        # check vector behavior
#        x, y = xf.focal_to_pixel([0, 100], [0, -100], 1, 1)
#        self.assertEqual(str([x, y]), "[array([ 1.,  6.]), array([ 1., -4.])]")

    def test_14_focal_to_pixel_degenerate(self):
        xf = Transformer()
        x, y, ccd_id = xf.focal_to_pixel(0, 0, 1,
                                              self.get_ccd_geometries(1))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([0]), array([0]), array([0])]")
        x, y, ccd_id = xf.focal_to_pixel(0, 0, 2,
                                              self.get_ccd_geometries(2))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([0]), array([0]), array([0])]")
        x, y, ccd_id = xf.focal_to_pixel(0, 0, 3,
                                              self.get_ccd_geometries(3))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([0]), array([0]), array([0])]")
        x, y, ccd_id = xf.focal_to_pixel(0, 0, 4,
                                              self.get_ccd_geometries(4))
        self.assertEqual(str([x, y, ccd_id]),
                         "[array([0]), array([0]), array([0])]")


if __name__ == '__main__':
    test_support.run_unittest(TransformerApiTestCase)
