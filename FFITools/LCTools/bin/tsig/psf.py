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
Loads and saves PSF data into and from a postgresql database.
"""

import os
import math
import logging
logger = logging.getLogger('tsig')

try:
    from scipy.io import loadmat
except ImportError:
    loadmat = None

import numpy as np

from tsig.util.db import Database

class SpreadFunctions(Database):
    """Loads PSF data from the database using queries"""
    CREDENTIALS_FILENAME = '~/.config/tsig/psf-dbinfo'
    CREDENTIALS_URL = 'http://tessellate.mit.edu/tsig/psf-dbinfo'

    @property
    def table(self):
        return self.dbinfo.get('dbtable', 'psf')

    def query(self, x, y, radius=3):
        """Make a query for nearby PSF data"""
        #
        # XXX This might need to change with the use of angles_x and angles_y instead of position point.
        #
        data = self._do_query("SELECT (camera_id, position, image, angles_x, angles_y) FROM " + self.table + " WHERE ST_DWithin(position, (%s, %s), %s)", [
            x, y, radius
        ])
        return (data['camera_id'], PSF(data['angles_x'], data['angles_y'], image))

    def add_to_db(self, **psf):
        """Add a single PSF to the database"""
        self._do_query("INSERT INTO " + self.table + " (type, temperature, angles_x, angles_y, position, image) VALUES (%s, %s, %s, %s, POINT(%s, %s), %s)", [
           str(psf['stellar_type'][:2]),
           int(psf['stellar_temp']),
           float(psf['field_angles'][0]),
           float(psf['field_angles'][1]),
           float(psf['field_position'][0]),
           float(psf['field_position'][1]),
           psf['PSFimage'].tolist(),
         ])  

    def add_file(self, filename, key='prf_bystellar'):
        """Insert all entries from a given matplot file"""
        if not loadmat:
            logger.error("SciPy is a required module to use PSF matlab loader")
            return

        if not os.path.isfile(filename):
            logger.error("File not found: %s" % filename)
            return

        try:
            mp = loadmat(filename)
        except Exception as err:
            logger.error("Error loading matplot file: %s" % str(err))
            return

        if key not in mp.keys():
            logger.error("Matplot file doesn't contain any PSFs")
            return

        for psf in mp[key][0]:
            psf = dict(zip(psf.dtype.names, psf))
            psf['stellar_temp'] = psf['stellar_temp'][0][0]
            psf['stellar_type'] = psf['stellar_type'][0][0][0]
            psf['field_angles'] = tuple(psf['field_angles'][0])
            psf['field_position'] = tuple(psf['field_position'][0])
            psf['PSFimage'] = psf['PSFimage']

            logger.debug("\n".join(
                ["Adding PSF %s: %s" % (name, str(psf[name])) for name in psf]))
            self.add_to_db(**psf)
        return True

    def create_table(self):
        """Convience function to create the known PSF database schema"""
        self._do_query("CREATE TABLE " + self.table + """(
  id          SERIAL     PRIMARY KEY,
  camera_id   INTEGER,
  type        CHAR(2),
  temperature INTEGER,
  angles_x    INTEGER,
  angles_y    INTEGER,
  position    POINT,
  image       REAL[][]
);""")


class PSF(object):
    """
    The origin for a PSF is lower left corner.
    """
    def __init__(self, field_angle_x, field_angle_y, grid=None):
        if grid is None:
            grid = np.zeros(shape=(1, 1), dtype=float)
        self.field_angle_x = field_angle_x
        self.field_angle_y = field_angle_y
        self.grid = grid

    def __repr__(self):
        """string representation of the PSF"""
        return self.get_label()

    def get_label(self):
        return "PSF fa_x=%s fa_y=%s w=%s h=%s sum=%.3f" % (
            self.field_angle_x, self.field_angle_y,
            self.grid.shape[0], self.grid.shape[1], self.grid.sum())

    def make_plot(self, ax=None, cmap='hot'):
        import matplotlib.pylab as plt
        if ax is None:
            plt.figure(self.get_label())
            ax = plt.subplot()
        else:
            plt.sca(ax)
        ax.imshow(self.grid, cmap=cmap, interpolation='none')
        ax.set_title(self.get_label())
        ax.invert_yaxis()
        return plt.gcf()

    def normalize(self):
        self.grid /= self.grid.sum()

    @staticmethod
    def _gradient(x, y, origin_x, origin_y, A=100.0):
        """
        Two-dimensional function with highest value at origin then diminishing
        in every direction from there.  Value drops off from maxval at the
        origin to zero, using a two-dimensional gaussian.

        origin is w,h
        """
        a = c = 0.001
        b = 0.0
        return A * math.exp(-(a * (x - origin_x) ** 2
                              - 2 * b * (x - origin_x) * (y - origin_y)
                              + c * (y - origin_y) ** 2))

    @staticmethod
    def create_random(field_angle_x=0.0, field_angle_y=0.0,
                      width=21, height=31):
        # shape is width, height
        return PSF(field_angle_x, field_angle_y,
                   np.random.random((width, height)))

    @staticmethod
    def create_gradient(field_angle_x=0.0, field_angle_y=0.0,
                        width=21, height=31, origin_x=0, origin_y=0):
        # shape is width, height
        # origin is width, height
        grid = np.zeros(shape=(height, width), dtype=float)
        for i in range(grid.shape[1]):
            for j in range(grid.shape[0]):
                grid[j, i] = PSF._gradient(
                    i, j, origin_x=origin_x, origin_y=origin_y)
        return PSF(field_angle_x, field_angle_y, grid)

    @staticmethod
    def interpolate(field_angle_x, field_angle_y, a, b, c, d):
        """
        Do a bilinear interpolation of 4 PSFs to obtain another PSF.

        x - x component of field angle
        y - y component of field angle
        a - PSF 0,0
        b - PSF 0,1
        c - PSF 1,0
        d - PSF 1,1

        Origin for x,y is lower left corner.
        """
        assert(a.grid.shape == b.grid.shape == c.grid.shape == d.grid.shape)
        psf = PSF(field_angle_x, field_angle_y, None)
        psf.grid = np.zeros(shape=a.grid.shape, dtype=float)
        x = field_angle_x
        y = field_angle_y
        delta_x = d.field_angle_x - a.field_angle_x
        delta_y = d.field_angle_y - a.field_angle_y
        wA = (d.field_angle_x - x) * (d.field_angle_y - y)
        wB = (x - c.field_angle_x) * (c.field_angle_y - y)
        wC = (b.field_angle_x - x) * (y - b.field_angle_y)
        wD = (x - a.field_angle_x) * (y - a.field_angle_y)
        for i in range(psf.grid.shape[1]):
            for j in range(psf.grid.shape[0]):
                v = wA * a.grid[j, i]
                v += wB * b.grid[j, i]
                v += wC * c.grid[j, i]
                v += wD * d.grid[j, i]
                v /= (delta_x * delta_y)
                psf.grid[j, i] = v
        return psf
