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
Access a TIC3 database and get stars based on position.
"""

import os
import time
import urllib
import numpy as np

from tsig.util.db import Database
from tsig.util import to_float
from .base import logger, Catalog

class TIC(Catalog, Database):
    NAME = "TIC3"
    CREDENTIALS_FILENAME = '~/.config/tsig/tic-dbinfo'
    CREDENTIALS_URL = 'http://tessellate.mit.edu/tsig/tic-dbinfo'

    TAG_RA = 'ra'
    TAG_DEC = 'dec'
    TAG_PMRA = 'pmra'
    TAG_PMDEC = 'pmdec'
    TAG_TMAG = 'tmag'
    TAG_TEFF = 'teff'
    COLUMNS = [TAG_RA, TAG_DEC, TAG_PMRA, TAG_PMDEC, TAG_TMAG, TAG_TEFF]

    def __init__(self, dbinfo_loc=None, dbtable='ticentries',
                 cachedir=None, cache_queries=True,
                 min_brightness=None, max_brightness=None, object_type=None,
                 row_limit=-1,
                 **dbinfo):
        Catalog.__init__(self)
        Database.__init__(self, dbinfo_loc=dbinfo_loc, dbtable=dbtable, **dbinfo)

        # query for everything, but put these limits on which objects we use
        self.min_brightness = to_float(min_brightness)
        self.max_brightness = to_float(max_brightness)
        self.object_type = object_type
        self.row_limit = row_limit

        # configure for caching
        if cachedir is None:
            cachedir = Catalog.get_cache_dirname()
        self.cachedir = cachedir
        Catalog.create_cache_dir(cachedir)
        self.cache_queries = cache_queries

    def query(self,
              ra=Catalog.DEFAULT_RA,
              dec=Catalog.DEFAULT_DEC,
              radius=Catalog.DEFAULT_RADIUS,
              dbindex=None):
        if dbindex == 'box':
            self.query_box_btree(ra, dec, 2 * radius, 2 * radius)
        elif dbindex == 'spacial_box':
            self.query_box(ra, dec, 2 * radius, 2 * radius)
        else:
            self.query_circle(ra, dec, radius)

    def query_circle(self,
                     ra=Catalog.DEFAULT_RA,
                     dec=Catalog.DEFAULT_DEC,
                     radius=Catalog.DEFAULT_RADIUS):
        """Do circle query using spacial index"""
        self.name = Catalog.create_query_name(self.NAME, ra, dec, radius=radius)
        sqlcmd = "\
SELECT %s FROM %s\
 WHERE spoint(radians(ra),radians(dec)) @ scircle '< (%sd,%sd), %sd >'" % (
     ','.join(self.COLUMNS), self.dbinfo['dbtable'], ra, dec, radius)
        self._query(sqlcmd)
        # remember the query
        self.query_ra = ra
        self.query_dec = dec
        self.query_radius = radius

    def query_box(self,
                  ra=Catalog.DEFAULT_RA,
                  dec=Catalog.DEFAULT_DEC,
                  width=Catalog.DEFAULT_WIDTH,
                  height=Catalog.DEFAULT_HEIGHT):
        """Do box query using spacial index"""
        self.name = Catalog.create_query_name(self.NAME, ra, dec,
                                              width=width, height=height)
        sqlcmd = "\
SELECT %s FROM %s\
 WHERE spoint(radians(ra),radians(dec)) @ sbox '( (%sd,%sd),(%sd,%sd) )'" % (
     ','.join(self.COLUMNS), self.dbinfo['dbtable'],
     ra - 0.5 * width, dec - 0.5 * height, ra + 0.5 * width, dec + 0.5 * height)
        self._query(sqlcmd)
        # remember the query
        self.query_ra = ra
        self.query_dec = dec
        self.query_width = width
        self.query_height = height

    def query_box_btree(self,
                        ra=Catalog.DEFAULT_RA,
                        dec=Catalog.DEFAULT_DEC,
                        width=Catalog.DEFAULT_WIDTH,
                        height=Catalog.DEFAULT_HEIGHT):
        """Do non-spacial box query"""
        self.name = Catalog.create_query_name(self.NAME, ra, dec,
                                              width=width, height=height)
        sqlcmd = "\
SELECT %s FROM %s\
 WHERE (ra between %s and %s) and (dec between %s and %s)" % (
     ','.join(self.COLUMNS), self.dbinfo['dbtable'],
     ra - 0.5 * width, ra + 0.5 * width, dec - 0.5 * height, dec + 0.5 * height)
        self._query(sqlcmd)
        # remember the query
        self.query_ra = ra
        self.query_dec = dec
        self.query_width = width
        self.query_height = height

    def _query(self, sqlcmd):
        """
        first try to load the catalog from local cache.  if that fails,
        make the query.
        """
        fn = self.name + '.npy'
        fullname = os.path.join(self.cachedir, fn)
        try:
            logger.debug("load stars from %s" % fullname)
            t = np.load(fullname)
        except IOError, e:
            logger.debug("load failed: %s" % e)
            self.check_dbinfo(self.dbinfo)
            logger.info("query %s: %s" % (self.NAME, sqlcmd))
            if self.row_limit is not None and self.row_limit > 0:
                sqlcmd += " LIMIT %s" % self.row_limit
            result = self._do_query(sqlcmd)
            dtype = [(self.COLUMNS[i], float) for i in range(len(self.COLUMNS))]
            t = np.array(result, dtype=dtype)
            if self.cache_queries:
                logger.debug("save query results to %s" % fullname)
                np.save(fullname, t)

        _ra = np.array(t[:][self.TAG_RA])
        _dec = np.array(t[:][self.TAG_DEC])
        _pmra = np.array(t[:][self.TAG_PMRA])
        _pmdec = np.array(t[:][self.TAG_PMDEC])
        _tmag = np.array(t[:][self.TAG_TMAG])
        _teff = np.array(t[:][self.TAG_TEFF])

        _pmra[np.isfinite(_pmra) == False] = 0.0
        _pmdec[np.isfinite(_pmdec) == False] = 0.0

        # FIXME: consider applying these constraints during query, not post
        ok = np.isfinite(_tmag)
        if self.min_brightness is not None:
            ok *= _tmag <= self.min_brightness
        if self.max_brightness is not None:
            ok *= _tmag >= self.max_brightness

        if len(_tmag):
            logger.debug("found %s objects; %s < Tmag < %s" %
                         (np.sum(ok), np.min(_tmag[ok]), np.max(_tmag[ok])))
        else:
            logger.debug("found 0 objects")

        self.ra = _ra[ok]
        self.dec = _dec[ok]
        self.pmra = _pmra[ok]
        self.pmdec = _pmdec[ok]
        self.tmag = _tmag[ok]
        self.teff = _teff[ok]
        self.epoch = 2000.0

    def query_by_id(self, tic_id, field_list=None):
        """Query the catalog by TIC identifier"""
        if isinstance(tic_id, int):
            tic_id = [tic_id]
        if field_list is None:
            field_list = '*'
        sqlcmd = "SELECT %s FROM %s WHERE id in (%s)" % (
            field_list, self.dbinfo['dbtable'], ','.join(tic_id))
        return self._do_query(sqlcmd)

    def query_by_loc(self, ra, dec, field_list=None):
        """Query the catalog by location ra,dec"""
        if field_list is None:
            field_list = '*'
        sqlcmd = "SELECT %s FROM %s WHERE ra=%s and dec=%s" % (
            field_list, self.dbinfo['dbtable'], ra, dec)
        return self._do_query(sqlcmd)
