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

__all__ = ('Spacecraft',)

import logging
import os

import numpy as np
import configobj

import tsig.camera
from tsig.util.configurable import ConfigurableObject, Configuration

from .transformer import Transformer


class Spacecraft(ConfigurableObject):
    """
    The spacecraft contains four cameras with pre-defined geometry.

    The spacecraft has a position and orientation.  The position (and velocity)
    matter only for differential velocity effects; they are not used in the
    coordinate transformations from starfield to CCD pixel.  The orientation
    is specified by the RA, Dec, and roll angle.

    A typical mission will capture data from a single sector for a two week
    period.  So for two weeks the RA, Dec, and roll angle will be constant
    as the spacecraft makes two orbits around earth.
    """
    SPACECRAFT_CFG = Configuration(
        os.path.join(os.path.dirname(__file__), 'spacecraft.cfg'))
    CAMERA_GEOMETRY_CFG = Configuration(
        os.path.join(os.path.dirname(__file__), 'camera-geometry.cfg'))
    CAMERA_FOCAL_CFG = Configuration(
        os.path.join(os.path.dirname(__file__), 'camera-focal.cfg'))
    CCD_GEOMETRY_CFG = Configuration(
        os.path.join(os.path.dirname(__file__), 'ccd-geometry.cfg'))
    CCD_PROPERTIES_CFG = Configuration(
        os.path.join(os.path.dirname(__file__), 'ccd-properties.cfg'))

    def __init__(self, **kwargs):
        """
        Create each camera using default geometry, override with any
        customizations from the configuration.
        """
        super(Spacecraft, self).__init__()

        self.data_set_id = kwargs.pop('data_set_id', None)
        self.geometry_model_id = kwargs.pop('geometry_model_id', None)
        self.start_tjd = kwargs.pop('start_tjd', None)

        # configure the cameras
        self.camera = []
        for i in range(1, 5):
            label = 'camera_%d' % i
            cam_dict = configobj.ConfigObj()
            cam_dict.merge(Spacecraft.SPACECRAFT_CFG.get(label))
            cam_dict.merge(Spacecraft.CAMERA_GEOMETRY_CFG.get(label))
            cam_dict.merge(Spacecraft.CAMERA_FOCAL_CFG.get(label))
            cam_dict.merge(Spacecraft.CCD_GEOMETRY_CFG.get(label))
            cam_dict.merge(Spacecraft.CCD_PROPERTIES_CFG.get(label))
            cam_dict['label'] = label
            cam_dict['number'] = i
            cam_dict.update(kwargs.get(label, {}))
            cam = tsig.camera.Camera(**cam_dict)
            self.camera.append(cam)

        # these are set when the camera is positioned and oriented
        self.ra = 0.0
        self.dec = 0.0
        self.roll = 0.0
        self.x = None # FIXME: need defaults for spacecraft position
        self.y = None
        self.z = None
        self.r = 0
        self.s = 0
        self.t = 0

    def __repr__(self):
        return "pointing=(%.2f,%.2f,%.2f) pos=(%.2f,%.2f,%.2f) vel=(%.2f,%.2f,%.2f)" % (
            self.ra, self.dec, self.roll,
            self.x, self.y, self.z, self.r, self.s, self.t)

    def set_position(self, x, y, z):
        """Position the spacecraft in orbit"""
        self.x = x
        self.y = y
        self.z = z

    def set_velocity(self, r, s, t):
        """Set the spacecraft velocity in space"""
        self.r = r
        self.s = s
        self.t = t

    def set_pointing(self, ra, dec, roll=0.0):
        """
        Aim the spacecraft at a point in the sky.

          ra   - Right ascension from the vernal equinox, in degrees
          dec  - Declination away from celestial equator, in degrees
          roll - The rotational angle of the spacecraft, in degrees

        """
        logging.info("point spacecraft to ra=%.2f dec=%.2f roll=%.2f" %
                     (ra, dec, roll))
        self.ra = ra
        self.dec = dec
        self.roll = roll
        xform = Transformer()
        for i in range(4):
            ra_cam, dec_cam, roll_cam = xform.spacecraft_to_camera_pointing(
                ra, dec, roll, self.camera[i].offset_angle)
            self.camera[i].point_camera(ra_cam, dec_cam, roll_cam)

    def get_cam_angle(self, cam_id):
        t = {1: -36.0, 2: -12.0, 3: 12.0, 4: 36.0}.get(cam_id)
        return t

    def get_cam_geometry(self, cam_id):
        return self.camera[cam_id - 1].get_geometry()

    def get_ccd_geometries(self, cam_id):
        return self.camera[cam_id - 1].get_ccd_geometries()

    def get_headers(self):
        return [
            ('SCRA', (self.ra, '[degree] spacecraft pointing right ascension')),
            ('SCDEC', (self.dec, '[degree] spacecraft pointing declination')),
            ('SCROLL', (self.roll, '[degree] spacecraft pointing roll angle')),
        ]
