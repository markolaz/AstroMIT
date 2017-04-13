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
A base class of test functionality for all possible tests.

We set up test data locations, caching locations and do clean up.
"""

import os
import sys
import shutil
import tempfile
import numpy

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '..'))

from unittest import TestCase
try:
    from test import test_support
except ImportError:
    from test import support as test_support


class TsigTestCase(TestCase):
    """
    Inherit from me to have access to test data, configurations and
    to have caching set up for testing.
    """
    @classmethod
    def setUpClass(cls):
        """Set up a test directory to store files"""
        # Set random number generator for predictable output
        numpy.random.seed(42)
        cls.test_dir = tempfile.mkdtemp(prefix='tsig')

    @classmethod
    def tearDownClass(cls):
        """Clean up created test directory"""
        if os.path.isdir(cls.test_dir):
	    shutil.rmtree(cls.test_dir)

    def setUp(self):
        """Create a copy of the test data and caching dir"""
        self.save_files = False
        self.save_cache = False

        # Set up test data directory
        self.data_dir = os.path.join(self.test_dir, 'data')
        shutil.copytree(os.path.join(DIR, 'data'), self.data_dir)

        # Set up caching
        self.cache_dir = os.path.join(self.test_dir, 'cache')
        os.makedirs(self.cache_dir)

    def tearDown(self):
	"""Remove copy of test data and caching dir"""
        if self.save_files:
            try:
                target = tempfile.mkdtemp(prefix='dumped-data-', dir=os.getcwd())
                os.rmdir(target)
                shutil.copytree(self.data_dir, target)
                sys.stdout.write("Kept data in: %s\n" % target)
            except Exception:
                sys.stdout.write("Failed to keep data.")

        if os.path.isdir(self.data_dir):
	    shutil.rmtree(self.data_dir)
        if os.path.isdir(self.cache_dir):
            shutil.rmtree(self.cache_dir)

    def assertDataFile(self, filename):
        """Return the filename of an existing data file"""
        filename = os.path.join(self.data_dir, filename)
        self.assertTrue(os.path.isfile(filename),
                        "Data file is missing: %s" % filename)
        return filename

    def assertCompareFile(self, filename):
        """Returns the filename and a a non-existant filename to populate
           so the original file can be compared.
        """
        filepath = self.assertDataFile(filename)
        newpath = "%s-output.%s" % tuple(filepath.rsplit('.', 1))
        self.assertFalse(os.path.isfile(newpath),
            "Output file exists before process to create it: %s" % newpath)
        return (filepath, newpath)

    def test_data(self):
        """Test data was created correctly"""
        self.assertTrue(self.data_dir)
        self.assertTrue(self.cache_dir)


