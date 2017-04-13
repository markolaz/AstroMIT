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
"""
# XXX - Turn a single long/lat into a full np.array of long lats

from .base import ImageEffect
import numpy as np

class BackgroundEffect(ImageEffect):
    def __init__(self, exptime=1800.0):
        super(BackgroundEffect, self).__init__()
        self.exptime = exptime

    def apply_to_image(self, image, calc):
	effective_area = 69.1
        pixelscale = 21.1
        pixel_solid_area = pixelscale ** 2
        bg = calc * effective_area * pixel_solid_area
        return image + bg * self.exptime


class CelestialBackground(BackgroundEffect):
    """
    Adds zodiacal light, treated as smooth using celestial long/lat
    """
    def __init__(self, longitude=0.0, latitude=0.0, exptime=1800.0):
        super(CelestialBackground, self).__init__(exptime=exptime)
        self.lon = longitude
        self.lat = latitude

    def apply_to_image(self, image):
        v_max = 23.345
        delta_v = 1.148
        b = np.abs(self.lat)
        v = v_max - delta_v * ((b - 90.0) / 90.0) ** 2
        assert ((b < 90).all())
        calc = 10 ** (-0.4 * (v - 22.8)) * (2.56e-3)
        return super(CelestialBackground, self).apply_to_image(image, calc)


class GalacticBackground(BackgroundEffect):
    """
    Adds unresolved stars, treated as smooth background using galactic long/lat
    """
    def __init__(self, longitude=0.0, latitude=0.0, exptime=1800.0):
        super(GalacticBackground, self).__init__(exptime=exptime)
        self.lon = longitude
        self.lat = latitude
        self.complicated = False

    def apply_to_image(self, image):

        # from Josh and Peter's memo
        flip = self.lon > 180.0
        if np.sum(flip):
            self.lon[flip] -= 360.0

        if self.complicated:
            L = np.abs(self.lon / 180.0)
            B = np.abs(self.lat / 90.0)
            a1 = 18.7
            a2 = 4.3
            a3 = 0.52
            a4 = 10.2
            a5 = 0.46
            a6 = -3.74
            I_surface_brightness = a1 + a2 * (1.0 - np.exp(-L / a3)) + a4 * (1.0 - np.exp(-B / a5)) + a6 * np.sqrt(
                L * B)
        else:
            a0 = 18.9733
            a1 = 8.833
            a2 = 4.007
            a3 = 0.805
            I_surface_brightness = a0 + a1 * (np.abs(self.lat) / 40.0) + a2 * (np.abs(self.lon) / 180.0) ** a3

        calc = 10 ** (-0.4 * I_surface_brightness) * 1.7e6
        return super(GalacticBackground, self).apply_to_image(image, calc)





