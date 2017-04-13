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

import logging
import numpy as np

from tsig.util.configurable import ConfigurableObject

logger = logging.getLogger(__name__)


class CCD(ConfigurableObject):
    """
    A CCD object used to emulate what a CCD chip will do when taking a photo.
    """

    def __init__(self, label, number=0, rows=2058, cols=2048, rotation=0.0,
                 x_0=0.0, y_0=0.0, pixel_size_x=0.015, pixel_size_y=0.015,
                 gain_a=1.0, gain_b=1.0, gain_c=1.0, gain_d=1.0,
                 saturation_limit=1000, # FIXME: need proper default
                 **kwargs):
        """
        Create a CCD object.

          label - A unique (within a single camera) label for this CCD.
          rows - number of pixel rows
          cols - number of pixel columns
          rotation - about the z axis
          x_0 - x origin offset, in mm
          y_0 - y origin offset, in mm
          pixel_size_x - size of each pixel, in mm
          pixel_size_y - size of each pixel, in mm
        """
        super(CCD, self).__init__()
        self.label = label
        self.number = number
        self.rows = int(rows)
        self.cols = int(cols)
        self.rotation = float(rotation)
        self.x_0 = float(x_0)
        self.y_0 = float(y_0)
        self.pixel_size_x = float(pixel_size_x)
        self.pixel_size_y = float(pixel_size_y)
        self.gain_a = float(gain_a)
        self.gain_b = float(gain_b)
        self.gain_c = float(gain_c)
        self.gain_d = float(gain_d)
        self.saturation_limit = int(saturation_limit)

    def __repr__(self):
        return self.label

    def get_geometry(self):
        return {
            'rows': self.rows,
            'cols': self.cols,
            'rotation': self.rotation,
            'x_0': self.x_0,
            'y_0': self.y_0,
            'pixel_size_x': self.pixel_size_x,
            'pixel_size_y': self.pixel_size_y
        }

    def get_headers(self):
        
        # HKBTMPC1 - ccd 1 board temp
        # HKBTMPC2 - ccd 2 board temp
        # HKBTMPC3 - ccd 3 board temp
        # HKBTMPC4 - ccd 4 board temp
        # HKDRVTMP - driver temperature
        # HKINTTMP - interface temperature
        # HKPTS01,02,03,04,05,06,07,08,09,10,11,12
        # HKALCUC1,C2,C3,C4 - ACLU sensor ccd N

        return [
            ('CCDLABEL', (self.label, '')),
            ('CCDNUM', (self.number, 'CCD number (1,2,3,4)')),
            ('CCDROWS', (self.rows, 'CCD rows')),
            ('CCDCOLS', (self.cols, 'CCD columns')),
            ('CCDXOFF', (self.x_0, '[mm] x offset from camera boresight')),
            ('CCDYOFF', (self.x_0, '[mm] y offset from camera boresight')),
            ('CCDPIXW', (self.pixel_size_x, '[mm] pixel width')),
            ('CCDPIXH', (self.pixel_size_y, '[mm] pixel height')),
            ('CCDROT', (self.rotation, '[degree] rotation around camera boresight')),
            ('CCDGAINA', (self.gain_a, '[electron/ADU] effective gain')),
            ('CCDGAINB', (self.gain_b, '[electron/ADU] effective gain')),
            ('CCDGAINC', (self.gain_c, '[electron/ADU] effective gain')),
            ('CCDGAIND', (self.gain_d, '[electron/ADU] effective gain')),
            ('CCDSAT', (self.saturation_limit, '[electron] saturation limit')),
        ]
