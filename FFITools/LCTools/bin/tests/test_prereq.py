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
Verify the environment and pre-requisites
"""

from unittest import TestCase

try:
    from test import test_support
except ImportError:
    from test import support as test_support

class PrereqTestCase(TestCase):

    def test_numpy(self):
        result = ""
        try:
            import numpy
            numpy.random.choice((0,1,2))
        except ImportError as e:
            result = str(e)
        except AttributeError as e:
            result = str(e)
        self.assertEqual(result, "")

    def test_astropy(self):
        result = ""
        try:
            import tempfile
            from astropy.io import fits
            filename = tempfile.TemporaryFile()
            hdulist = fits.HDUList()
            hdulist.writeto(filename, overwrite=True)
        except TypeError as e:
            result = str(e)
        self.assertEqual(result, "")


if __name__ == '__main__':
    test_support.run_unittest(PrereqTestCase)
