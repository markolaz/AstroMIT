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
Loads and saves PRNU data into and from a postgresql database.
"""

import os
import re
import math
import logging
logger = logging.getLogger('tsig')

import numpy as np
from astropy.io import fits

from tsig.util.db import Database

class PhotoResponse(Database):
    """Loads PSF data from the database using queries"""
    CREDENTIALS_FILENAME = '~/.config/tsig/prnu-dbinfo'
    CREDENTIALS_URL = 'http://tessellate.mit.edu/tsig/prnu-dbinfo'

    # This matches the current prnu file specification
    # d - Detector level, always d for prnu. If c, then camera level and not prnu.
    # board - Board version, usually 7
    # cam - Camera Serial, e.g. 05
    # series - If more than one set was taken for the same settings, this increases by one
    # wave - broadband or flood waveband number
    # temp - Temperature of the exposure
    # inter - o32 -> 2d 32nd order legendre polynomial fit to the illumination function
    FN_MATCH = r'd(?P<board>\d)(?P<cam>\d{2})(?P<series>\d{2})_(flood)?(?P<wave>(bband|\d+))_m(?P<temp>\d+)_prnu_(?P<inter>\w+)\.fits'

    @property
    def table(self):
        return self.dbinfo.get('dbtable', 'prnu')

    def query(self, camera, wavelength='bband', temperature=70):
        """Make a query for PRNU data"""
        data = self._do_query("SELECT fits FROM " + self.table + \
            " WHERE camera_id=%s AND wavelength=%s AND temperature=%d", [
            camera, wavelength, temperature,
        ])
        return [np.array(datum['fits']) for datum in data]

    def add_to_db(self, **prnu):
        """Add a single PRNU to the database"""
        self._do_query("INSERT INTO " + self.table + " (ccd, board, camera, series, wavelength, temperature, interpolation, fits) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)", [
           int(prnu['ccd']),
           int(prnu['board']),
           int(prnu['cam']),
           int(prnu['series']),
           int(prnu['wave']),
           int(prnu['temp']),
           prnu['inter']
           prnu['Fits'].tolist(),
         ])  

    def add_file(self, filename):
        """Insert all entries from a given fits file"""
        if not os.path.isfile(filename):
            logger.error("File not found: %s" % filename)
            return

        name = os.path.basename(filename)
        match = re.search(self.FN_MATCH, name)
        if not match:
            logger.error("Filename not in PRNU filename format: %s" % name)
            return

        data = match.groupdict()
        if data['wave'] == 'bband':
            data['wave'] = 0

        try:
            logging.debug("Loading fits file: %s" % filename)
            hdulist = fits.open(filename, memmap=True)
            logging.debug("DONE")
        except Exception as err:
            logger.error("Error loading fits file: %s" % str(err))
            return

        img = hdulist[0].data
        height_a, width_a = img.shape
        width_b, height_b = (4272, 4156)
        if width_a != width_b or height_b != height_b:
            logging.error("PRNU fits file isn't the right size. Found "
                "%dx%d expected %dx%d" % (width_a, height_a, width_b, height_b))
            return

        left = (0, 44)
        right = (4229, 4272)
        top = (-1, -1)
        bottom = (5156, 5156)
        mid_x = (2093, 2180)
        # Original coordinates give larger mid_y band (horizontal bar)
        #mid_y = (2049, 2108)
        # But we need a slightly larger range.
        mid_y = (2059, 2098)

        for i, (x1, x2) in enumerate([(left[1], mid_x[0]), (mid_x[1], right[0])]):
            for j, (y1, y2) in enumerate([(top[1], mid_y[0]), (mid_y[1], bottom[0])]):
                logging.debug("Taking section %dx%d-%dx%d" % (x1+1, y1+1, x2, y2))
                # The CCD number is re-mapped from a1,b1,a2,b2 to
                # anti-clockwise from top-right, See CCD reference numbers.
                data['ccd'] = [2,1,3,4][i + (j * 2)]
                logging.debug("Prepared CCD %(ccd)s for Camera %(cam)s" % data)
                # The fits image we save is just the science bits for one ccd
                data['Fits'] = img[y1+1:y2,x1+1:x2]
                logging.debug("New shape: %dx%d" % data['Fits'].shape)
                self.add_to_db(**data)
        return True

    def create_table(self):
        """Convience function to create the known PSF database schema"""
        self._do_query("CREATE TABLE " + self.table + """(
  id          SERIAL     PRIMARY KEY,
  ccd         INTEGER,
  camera      INTEGER,
  temperature INTEGER,
  wavelength  INTEGER,
  series      INTEGER,
  board       INTEGER,
  interpolation CHAR(12),
  fits        REAL[][]
);""")


