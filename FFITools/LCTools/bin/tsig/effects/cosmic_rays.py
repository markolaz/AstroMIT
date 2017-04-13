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
This effect requires that tsig_cosmical is installed (C code module)
"""

from .base import ImageEffect
import numpy as np

class CosmicRays(ImageEffect):
    """
    Noise effect which deposits a track of charge across the pixels depending
    on the direction. Modified lamdel distribution.
    """
    def __init__(self, exptime=1800.0):
        """
        exptime - Exposure time on the CCD
        """
        super(CosmicRays, self).__init__()
        self.exptime = exptime
        self.rate = 5.0
        self.gradient = False
        self.diffusion = False
        self.buffer = 100

    def apply_to_image(self, image):
        try:
            from tsig_cosmical import cosmical
        except ImportError:
            raise ImportError("The tsig_cosmical extension must be installed to use the Cosmic Ray Effect.")

        size = max(image.shape)
        smallexptime = self.exptime * 1.5 
        bigexptime = self.exptime * 1.5 
        # if a gradient is set, allow the exposure times to be different
        if self.gradient:
            smallexptime = exptime
            bigexptime = exptime * 2 

        # add margin around image, because Al's code doesn't start cosmic ray
        # images off the screen
        bufferedsize = size + 2 * self.buffer

        # set the diffusion flag to be 0 if False, 1 if True (as integer type)
        intdiffusion = 0  # np.int(self.diffusion)

        # call the fancy cosmic ray code
        bg = cosmical(self.rate, smallexptime, bigexptime,
                      bufferedsize, bufferedsize, intdiffusion)

        # if we need to diffuse the image, use Al's kernal (from the fancy code)
        if self.diffusion:
            import scipy.signal
            kernal = np.array([[0.0034, 0.0516, 0.0034],
                               [0.0516, 0.7798, 0.0516],
                               [0.0034, 0.0516, 0.0034]])
            bg = scipy.signal.convolve2d(bg, kernal, mode='same')

        # Return the resized image on top of the fits file image
        buf = self.buffer
        img = bg[buf:-buf, buf:-buf] if buf > 0 else bg[:, :]
        return image + img

