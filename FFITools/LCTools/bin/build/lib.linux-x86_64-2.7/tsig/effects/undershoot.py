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
import numpy as np

class Undershoot(ImageEffect):
    """
    Apply an Undershoot effect
    """
    def __init__(self, amount=0.0013):
	self.amount = amount

    def apply_to_image(self, image):
        """
        Add undershoot to an image
        """
	kernel = np.array([1.0, -self.amount])
	def convolve_row(row):
	    return np.convolve(row, kernel, mode='same')
        return image + np.apply_along_axis(convolve_row, 1, image)

