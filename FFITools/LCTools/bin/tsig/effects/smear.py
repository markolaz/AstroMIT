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
A smear effect which can be applied to a CCD() object to simulate
how new pixels are exposed as the pixels as read in row by row.
"""

from .base import ImageEffect
import numpy as np

class Smear(ImageEffect):
    """
    Create a smear effect which can be applied to a CCD() object.
    """
    def __init__(self, readout_time=0.02, single_read=2.0):
        """
          readout_time - The amount of time it takes read a single row
          single_read - Amount of time the ccd read the image
        """
        self.single_read = float(single_read) # seconds
        self.readout_time = float(readout_time) # seconds

    def apply_to_image(self, image):
        """
        Apply a smear to an image.
        """
        # FIXME: default to single_read from EXPTIME in image headers.  if that
        # does not exist or if we have a value, then use our value
        ones = np.ones(image.shape)
        # Select each column, creating a mean of the column
        mean = np.mean(image, 0).reshape(1, image.shape[1]) * ones
        image += mean * self.readout_time / self.single_read
        return image

    def get_headers(self):
        return [
            ('SMEARROT', (self.readout_time, '[s] readout time')),
            ('SMEARSRT', (self.single_read, '[s] single read time')),
        ]
