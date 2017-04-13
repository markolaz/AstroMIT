#
# Copyright (C) 2015 - Zach Berta-Thompson <zkbt@mit.edu> (MIT License)
#               2017 - Massachusetts Institute of Technology (MIT)
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
"""Transformations from one reference frame to another.

cam_id - camera identifier, one of 1, 2, 3, or 4
ccd_id - CCD identifier, one of 1, 2, 3, or 4
"""

import astropy
import math
import numpy as np


class Transformer(object):
    """
    This object knows how to transform from one reference frame to another.
    """
    def spacecraft_to_camera_pointing(self, ra_sc, dec_sc, roll_sc, cam_angle):
        """
        Get the (ra,dec,roll) of a camera given the spacecraft (ra,dec,roll).
        """
        # get the tranformation matrix from equatorial to camera
        m1 = self.get_eq_sc_matrix(ra_sc, dec_sc, roll_sc)
        m2 = self.get_sc_cam_matrix(cam_angle)
        m = np.dot(m2, m1)
        # unit vectors in the camera reference frame
        x_cam = m[0]
        y_cam = m[1]
        z_cam = m[2]
        # unit vector of z axis in the equatorial
        z_eq = np.array([0, 0, 1])
        # from those get the ra and dec
        ra_cam = math.atan2(z_cam[1], z_cam[0]) * 180.0 / math.pi
        dec_cam = 90.0 - math.acos(z_cam[2]) * 180.0 / math.pi
        # get the east and north unit vectors to get roll
        e_vec = np.cross(z_eq, z_cam) / np.linalg.norm(np.cross(z_eq, z_cam))
        n_vec = np.cross(z_cam, e_vec)
        x = np.dot(x_cam, e_vec)
        y = np.dot(x_cam, n_vec)
        roll_cam = math.atan2(y, x) * 180.0 / math.pi
        # FIXME: deal with case where z_eq x z_cam is zero
        return ra_cam, dec_cam, roll_cam

    def celestial_to_focal(self, ra, dec, ra_sc, dec_sc, roll_sc, cam_geometry):
        """
        Transform celestial coordinates (ra,dec) to focal plane (x,y)
        on indicated camera with spacecraft oriented at (ra_sc,dec_sc,roll_sc)
        """
        # transform ra,dec to camera ra,dec
        ra_cam, dec_cam, fa_x, fa_y = self.celestial_to_camera(
            ra, dec, ra_sc, dec_sc, roll_sc, cam_geometry['angle'])
        # transform ra,dec to focal coordinates
        x, y = self.camera_to_focal(ra_cam, dec_cam, cam_geometry)
        return x, y, fa_x, fa_y

    def celestial_to_pixel(self, ra, dec, ra_sc, dec_sc, roll_sc,
                           cam_geometry, ccd_geometries, subpixels_per_pixel=1):
        """
        Transform celestial coordinates (ra,dec) to CCD pixel (col,row)
        on indicated camera with spacecraft oriented at (ra_sc,dec_sc,roll_sc)
        """
        # transform ra,dec to camera ra,dec
        ra_cam, dec_cam, fa_x, fa_y = self.celestial_to_camera(
            ra, dec, ra_sc, dec_sc, roll_sc, cam_geometry['angle'])
        # transform ra,dec to focal coordinates
        x, y = self.camera_to_focal(ra_cam, dec_cam, cam_geometry)
        # transform ra,dec to pixels
        x_n, y_n, ccd = self.focal_to_pixel(x, y,
                                            cam_geometry['number'],
                                            ccd_geometries,
                                            subpixels_per_pixel)
        return x_n, y_n, ccd, fa_x, fa_y

    def celestial_to_camera(self, ra, dec, ra_sc, dec_sc, roll_sc, cam_angle):
        """
        Transform celestial (ra,dec) to camera (ra,dec) given spacecraft
        pointing (ra_sc,dec_sc,roll_sc).  Also return field angle in two
        components fa_x and fa_y.
        """
        # ensure that the inputs are arrays
        if isinstance(ra, (type(None), int, float)):
            ra = [ra]
        if isinstance(dec, (type(None), int, float)):
            dec = [dec]
        # allocate ra,dec arrays for the results
        ra_cam = np.array(ra, dtype=float)
        dec_cam = np.array(dec, dtype=float)
        fa_x = np.zeros(shape=ra_cam.shape)
        fa_y = np.zeros(shape=ra_cam.shape)
        # generate the transformation matrix
        m1 = self.get_eq_sc_matrix(ra_sc, dec_sc, roll_sc)
        m2 = self.get_sc_cam_matrix(cam_angle)
        m = np.dot(m2, m1)
        # populate the ra,dec arrays by applying the transform to each element
        for i in range(ra_cam.size):
            ra_r = ra_cam[i] * math.pi / 180.0
            dec_r = dec_cam[i] * math.pi / 180.0
            v = np.array([math.cos(ra_r) * math.cos(dec_r),
                          math.sin(ra_r) * math.cos(dec_r),
                          math.sin(dec_r)])
            c = np.dot(m, v)
            ra_cam[i] = math.atan2(c[1], c[0]) * 180.0 / math.pi
            try:
                dec_cam[i] = math.asin(c[2]) * 180.0 / math.pi
            except ValueError:
                dec_cam[i] = 90.0
            if c[2]:
                fa_x[i] = math.atan(c[0] / c[2]) * 180.0 / math.pi
                fa_y[i] = math.atan(c[1] / c[2]) * 180.0 / math.pi
            else:
                fa_x[i] = 90.0
                fa_y[i] = 90.0
        return ra_cam, dec_cam, fa_x, fa_y

    def camera_to_focal(self, ra, dec, cam_geometry):
        """
        Transform a camera celestial (ra,dec) to (x,y) on the focal plane.
        Return value is in mm.

        ra - right ascension in the camera refrence frame
        dec - declination in the camera reference frame
        camera - camera number 1,2,3, or 4
        ccd - CCD number 1,2,3, or 4
        """
        lon = np.array(ra)
        colat = 90.0 - np.array(dec)
        t = np.tan(colat * math.pi / 180.0)
        t2 = t * t
        geom = cam_geometry
        # FIXME: include higher-order terms
        r = geom['focal_length'] * t * (geom['a_0'] + geom['a_2'] * t2 + geom['a_4'] * t2 ** 2 + geom['a_6'] * t2 ** 3 + geom['a_8'] * t2 ** 4)
        x = -r * np.cos(lon * math.pi / 180.0)
        y = -r * np.sin(lon * math.pi / 180.0)
        return x, y

    def focal_to_pixel(self, x, y, cam_number, ccd_geometries,
                       subpixels_per_pixel=1):
        """
        Transform (x,y) on the focal plane to pixel (col,row) on the CCD.

        x, y - coordinates in the focal plane
        camera - camera number 1,2,3, or 4

        return:
        x_n, y_n - pixel coordinates on the CCD
        ccd - CCD number 1,2,3, or 4

        geom - geometric properties for a specific camera and ccd
        geom['rotation'] - CCD rotation about the camera boresight
        geom['x_0'] - bottom left x coordinate of imaging area, in mm
        geom['y_0'] - bottom left y coordinate of imaging area, in mm
        geom['pixel_size_x'] - pixel dimension, in mm
        geom['pixel_size_y'] - pixel dimension, in mm

        x_fpp, y_fpp - in the fpp frame, the +x axis runs from the origin at
          the center of the array of 4 CCDs between CCDs 1 and 4, and the +y
          axis runs from the origin between CCDs 1 and 2.

        x_fpr, y_fpr - frame in which x,y axes are parallel to the CCD rows
          and columns, but rotated relative to the fpp frame
        """
        from math import cos, sin
        if isinstance(x, (type(None), int, float)):
            x = [x]
        if isinstance(y, (type(None), int, float)):
            y = [y]
        x_fpp = np.array(x)
        y_fpp = np.array(y)
        if cam_number in [1, 2]:
            x_fpp *= -1
            y_fpp *= -1
        x_n = np.zeros(shape=x_fpp.shape, dtype=int)
        y_n = np.zeros(shape=y_fpp.shape, dtype=int)
        ccd_id = np.zeros(shape=x_fpp.shape, dtype=int)
        for i in range(x_fpp.size):
            if x_fpp[i] > 0 and y_fpp[i] > 0:
                ccd_id[i] = 1
            elif x_fpp[i] < 0 and y_fpp[i] > 0:
                ccd_id[i] = 2
            elif x_fpp[i] < 0 and y_fpp[i] < 0:
                ccd_id[i] = 3
            elif x_fpp[i] > 0 and y_fpp[i] < 0:
                ccd_id[i] = 4
            if ccd_id[i]:
                g = ccd_geometries[ccd_id[i]]
                r = g['rotation'] * math.pi / 180.0
                x_delta = x_fpp[i] - g['x_0']
                y_delta = y_fpp[i] - g['y_0']
                x_fpr = x_delta * cos(r) + y_delta * sin(r)
                x_fpr *= subpixels_per_pixel
                y_fpr = y_delta * cos(r) - x_delta * sin(r)
                y_fpr *= subpixels_per_pixel
                x_n[i] = 1 + np.rint(x_fpr / g['pixel_size_x'])
                y_n[i] = 1 + np.rint(y_fpr / g['pixel_size_y'])
        return x_n, y_n, ccd_id

    def get_eq_sc_matrix(self, ra, dec, roll):
        """Matrix to transform from celestial equatorial to spacecraft"""
        from math import cos, sin
        ra_r = ra * np.pi / 180.0
        dec_r = dec * np.pi / 180.0
        roll_r = roll * np.pi / 180.0
        return np.array(
            [[-cos(roll_r) * sin(dec_r) * cos(ra_r) + sin(roll_r) * sin(ra_r),
              -cos(roll_r) * sin(dec_r) * sin(ra_r) - sin(roll_r) * cos(ra_r),
              cos(roll_r) * cos(dec_r)],
             [sin(roll_r) * sin(dec_r) * cos(ra_r) + cos(roll_r) * sin(ra_r),
              sin(roll_r) * sin(dec_r) * sin(ra_r) - cos(roll_r) * cos(ra_r),
              -sin(roll_r) * cos(dec_r)],
             [cos(dec_r) * cos(ra_r),
              cos(dec_r) * sin(ra_r),
              sin(dec_r)]])

    def get_sc_cam_matrix(self, cam_angle):
        """Matrix to transform from spacecraft to camera coordinates"""
        t = cam_angle * math.pi / 180.0
        return np.array([[0, -1, 0],
                         [math.cos(t), 0, -math.sin(t)],
                         [math.sin(t), 0, math.cos(t)]])
