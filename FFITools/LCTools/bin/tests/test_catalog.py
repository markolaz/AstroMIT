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
Test basic catalog functions using the TestGrid catalog.
"""

from unittest.case import SkipTest
from base import TsigTestCase, test_support

from tsig.catalog import *

class CatalogTestCase(TsigTestCase):
    def test_00_test_pattern(self):
        """Check the default test pattern catalog"""
        c = TestGrid()
        c = TestGrid({'TestGrid': {'size': 100, 'spacing': 10}})
        c = TestGrid({'TestGrid': {'center': (0.5, 1.5)}})
        c = TestGrid({'TestGrid': {'mag_range': (1, 5)}})
        c = TestGrid({'TestGrid': {'max_nudge': 1.2}})
        c = TestGrid({'TestGrid': {'max_motion': 1.2}})
        c = TestGrid({'TestGrid': {'randomize_magnitudes': True}})

    def test_02_constant_lightcurves(self):
        """Default should have all constant curves"""
        c = TestGrid({'TestGrid': {'size': 10, 'spacing': 5}})
        c.query()
        c.add_lightcurves()
        for x in c.lightcurve_codes:
            assert(x.startswith('ConstantLightCurve(offset'))

    def test_03_test_pattern_static(self):
        """Verify static data for test pattern"""
        c = TestGrid({'TestGrid': {'size': 10, 'spacing': 1}})
        c.stars_static()

    def test_04_test_pattern_snapshot(self):
        """Verify snapshot for test pattern"""
        c = TestGrid({'TestGrid': {'size': 10, 'spacing': 1}})
#        c.stars_snapshot(bjd=1) # FIXME: what value for bjd?
#        c.stars_snapshot(epoch=2001) # FIXME: what value for epoch?
#        c.stars_snapshot(exptime=1.0) # FIXME: what value for exptime?
#        c.stars_snapshot(roll=1.0) # FIXME: what value for roll?
#        c.stars_snapshot(ra_0=1.0, dec_0=1.0) # FIXME: what values?

    def test_05_test_pattern_plot(self):
        try:
            import matplotlib
            # TkAgg will show plots as windows, Agg will allow png comparisons.
            matplotlib.use('Agg')
        except ImportError:
            raise SkipTest("Matplotlib not installed.")

        from matplotlib.testing.compare import compare_images
        from matplotlib import pyplot as plt

        for c in [TestGrid(),
                  TestGrid({'TestGrid': {'size': 100, 'spacing': 10}}),
                  TestGrid({'TestGrid': {'size': 1000, 'spacing': 50}}),
                  TestGrid({'TestGrid': {'center': (0.5, 1.5)}})]:
            c.query()
            c.add_lightcurves()
            figure = c.make_plot()
            figure.canvas.draw_idle()
#            (f1, f2) = self.assertCompareFile('catalog/%s.png' % type(c).__name__)
#            figure.savefig('catalog_output.png')
            plt.close('all')

#            result = compare_images(f1, f2, 0.01)
#            if result:
#                self.save_files = True
#            self.assertFalse(result, str(result))


if __name__ == '__main__':
    test_support.run_unittest(CatalogTestCase)
