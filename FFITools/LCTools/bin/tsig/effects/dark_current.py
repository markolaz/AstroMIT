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

from .base import ImageEffect
import math
import numpy as np

class DarkCurrent(ImageEffect):
    """
    Apply dark noise to an object.

    For each port, the average dark current is measured in electrons per
    pixel per second.

    Dark current is modeled as an offset plus a gaussian perturbation to that
    offset where the mean of the gaussian is one over the square root of the
    offset.

    For a 2 second exposure, a typical dark current will be around 5 electrons
    per pixel.

    reference: Deb Woods 23feb2017
    """
    def __init__(self, offset=0.08):
        super(DarkCurrent, self).__init__()
        self.offset = 0.08 # electrons per pixel per second
        self.mean = 1.0 / math.sqrt(self.offset)

    def apply_to_image(self, image):
        image += np.random.normal(0, self.mean, len(image)) + self.offset
        return image

