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
Add different noise effects to an image slice.
"""

import numpy as np

from .base import ImageEffect

class ReadoutNoise(ImageEffect):
    """Apply Read Noise to an object.

    Read noise is modeled by apply a gaussian perturbation to each pixel.
    The mean value of the gaussian is different for each of the 4 ports in
    each of the CCDs.  Values are provided by calibration tables, for example
    
    "SN06 basic calibration at -70C"
    https://tess-web.mit.edu/ground/PNY-1/d70601/summaries/basicCal-70C

    reference: Deb Woods 23feb2017
    """
    def __init__(self, mean=100.0):
        super(ReadoutNoise, self).__init__()
        self.mean = mean

    def apply_to_image(self, image):
        # This should take account of the read size and take the
        # mean of the image value from the source data correctly.
        image += np.random.normal(0, self.mean, len(image))
        return image

    def get_headers(self):
        return [
            ('READNOIS',
             (self.mean, '[e-] read noise (per individual read)')),
        ]


class ShotNoise(ImageEffect):
    """Also known as Photon Noise"""
    def __init__(self):
        super(ShotNoise, self).__init__()

    def apply_to_image(self, image):
        ok = image > 0
        noise = np.zeros_like(image)
        noise[ok] = np.sqrt(image[ok])
        image += noise * np.random.randn(len(image))
        return image

    def get_headers(self):
        # FIXME: what parameter and units for shot noise header?
        return [
            ('SHOTNOIS', (0, '[e-] shot noise')),
        ]


class FixedPatternNoise(ImageEffect):
    def __init__(self):
        super(FixedPatternNoise, self).__init__()

    def apply_to_image(self, image):
        raise NotImplementedError("Fixed Pattern Noise Effect not written yet.")


class BlackSky(ImageEffect):
    def __init__(self):
        super(BlackSky, self).__init__()

    def apply_to_image(self, image):
        raise NotImplementedError("Black Sky Effect not written yet.")


