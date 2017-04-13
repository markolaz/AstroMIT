#!/usr/bin/env python
#
# Copyright (C) 2016 - Massachusetts Institute of Technology (MIT)
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
Test the TIC.
"""

import os
import sys

from unittest.case import SkipTest
from base import TsigTestCase, test_support

from tsig.catalog import *

class TICTestCase(TsigTestCase):

    RA = 30.0
    DEC = 50.0
    RADIUS = 0.2

    def setUp(self):
        super(TICTestCase, self).setUp()
        try:
            dbinfo = TIC.read_dbinfo()
            self._cfg = {'TIC': dbinfo}
        except Exception as e:
            raise SkipTest(e.message)

    def test_01_default_query(self):
        """Check default TIC"""
        c = TIC(self._cfg)
        c.query()

    def test_02_btree_indexing(self):
        """query with btree index"""
        import time
        self._cfg['TIC']['cache_queries'] = False
        c = TIC(self._cfg)
        start = time.time()
        c.query(ra=self.RA, dec=self.DEC, radius=self.RADIUS, dbindex='btree')
        elapsed = time.time() - start
        print "btree: %d records in %.3fs" % (len(c.tmag), elapsed)

    def test_02_box_indexing(self):
        """query with spacial index box"""
        import time
        self._cfg['TIC']['cache_queries'] = False
        c = TIC(self._cfg)
        start = time.time()
        c.query(ra=self.RA, dec=self.DEC, radius=self.RADIUS, dbindex='spacial_box')
        elapsed = time.time() - start
        print "box: %d records in %.3fs" % (len(c.tmag), elapsed)

    def test_02_circle_indexing(self):
        """query with spacial index circle"""
        import time
        self._cfg['TIC']['cache_queries'] = False
        c = TIC(self._cfg)
        start = time.time()
        c.query(ra=self.RA, dec=self.DEC, radius=self.RADIUS, dbindex='spacial_circle')
        elapsed = time.time() - start
        print "circle: %d records in %.3fs" % (len(c.tmag), elapsed)

if __name__ == '__main__':
    test_support.run_unittest(TICTestCase)
