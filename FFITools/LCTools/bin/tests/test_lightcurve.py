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
Test light curve base classes.
"""

import os
import sys
import numpy.random
import numpy as np

from unittest.case import SkipTest
from base import TsigTestCase, test_support

from tsig.lightcurve import *

class LightCurveTestCase(TsigTestCase):
    """Unit tests for the light curve classes."""
    def test_00_load(self):
        """Can the light curve object be instantiated"""
        self.assertTrue(isinstance(LightCurve(), dict))

    def test_01_constant(self):
        """Ask for a constant light curve"""
        ConstantLightCurve()

    def test_02_trapezoid(self):
        """Ask for a trapezoid light curve with default options"""
        TrapezoidLightCurve()

    def test_03_sinusoid(self):
        """Ask for a sinusoid light curve with default options"""
        SinusoidLightCurve()

    def test_04_mcquillan(self):
        """Ask for a McQuillan light curve with default options"""
        McQuillanLightCurve()

    def test_05_keplertce(self):
        """Ask for a KeplerTCE light curve with default options"""
        KeplerTCELightCurve()

    def test_10_random(self):
        """Ask for a random light curve"""
        lightcurve.create_random_lightcurve()

    def test_11_extreme(self):
        """Ask for an extreme light curve"""
        lightcurve.create_cartoon_lightcurve()

    def test_12_plot(self):
        try:
            import matplotlib
            # TkAgg will show plots as windows, Agg will allow png comparisons.
            matplotlib.use('Agg')
        except ImportError:
            raise SkipTest("Matplotlib not installed.")

        from matplotlib.testing.compare import compare_images
        from matplotlib import pyplot as plt

        for cls in [ConstantLightCurve,
                    TrapezoidLightCurve,
                    SinusoidLightCurve,
                    McQuillanLightCurve,
                    KeplerTCELightCurve]:
            obj = cls()
            figure = obj.make_plot()
            figure.canvas.draw_idle()
            (f1, f2) = self.assertCompareFile('lightcurves/%s.png' % cls.__name__)
            figure.savefig(f2)
            plt.close('all')

            result = compare_images(f1, f2, 0.01)
            if result:
                self.save_files = True
            self.assertFalse(result, str(result))

    def test_12_integration(self):
        """Test the integration"""

        # value of t does not matter
        c = ConstantLightCurve()
        x = c.integrated(0)
        self.assertEqual(x, 0)
        c = ConstantLightCurve()
        x = c.integrated(1)
        self.assertEqual(x, [0])
        c = ConstantLightCurve()
        x = c.integrated(10)
        self.assertEqual(x, [0])

        # check type safety of the arguments
        c = ConstantLightCurve(offset=1.0)
        x = c.integrated(np.array([1,2,3]))
        self.assertEqual(str(x), '[ 1.  1.  1.]')
        c = ConstantLightCurve(offset=1.0)
        x = c.integrated([1,2,3])
        self.assertEqual(str(x), '[ 1.  1.  1.]')
        c = ConstantLightCurve(offset=1.0)
        x = c.integrated((1,2,3))
        self.assertEqual(str(x), '[ 1.  1.  1.]')

        # check the constant offset
        c = ConstantLightCurve(offset=10.0)
        x = c.integrated(1)
        self.assertEqual(x, [10.0])
        c = ConstantLightCurve(offset=10.0)
        x = c.integrated(np.array([1,2,3]))
        self.assertEqual(str(x), '[ 10.  10.  10.]')

        # exposure time
        c = ConstantLightCurve(offset=1.0)
        x = c.integrated(1, 2.0 / 60.0 / 24.0)
        self.assertEqual(x, [1.])
        c = ConstantLightCurve(offset=1.0)
        x = c.integrated(1, 30.0 / 60.0 / 24.0)
        self.assertEqual(x, [1.])

        # check a sinusoid curve
        c = SinusoidLightCurve()
        x = c.integrated(1)
        self.assertEqual(str(x), '[ 0.09092291]')
        c = SinusoidLightCurve()
        x = c.integrated(np.array([1,2,3]))
        self.assertEqual(str(x), '[ 0.09092291 -0.07567495 -0.02794011]')

        # exposure time
        c = SinusoidLightCurve()
        x = c.integrated(1, 2.0 / 60.0 / 60.0 / 24.0)
        self.assertEqual(str(x), '[ 0.09092974]')
        c = SinusoidLightCurve()
        x = c.integrated(1, 2.0 / 60.0 / 24.0)
        self.assertEqual(str(x), '[ 0.09092971]')
        c = SinusoidLightCurve()
        x = c.integrated(1, 30.0 / 60.0 / 24.0)
        self.assertEqual(str(x), '[ 0.09092291]')

        # resolution has no effect on sinusoid
        c = SinusoidLightCurve()
        x = c.integrated(1, 30.0 / 60.0 / 60.0 / 24.0, 10)
        self.assertEqual(str(x), '[ 0.09092974]')
        c = SinusoidLightCurve()
        x = c.integrated(1, 2.0 / 60.0 / 60.0 / 24.0, 1000)
        self.assertEqual(str(x), '[ 0.09092974]')
        c = SinusoidLightCurve()
        x = c.integrated(1, 2.0 / 60.0 / 60.0 / 24.0, 10000)
        self.assertEqual(str(x), '[ 0.09092974]')

class RandomLightCurveTestCase(TsigTestCase):

    def test_00_randomness(self):
        rlc = RandomLightCurve()
        a = [rlc.get_flux(t*30) for t in xrange(10)]
        b = [rlc.get_flux(t*30) for t in xrange(10)]
        self.assertNotEqual(a, b)

    def test_01_randomseed(self):
        random.seed(100)
        rlc = RandomLightCurve()
        a = [rlc.get_flux(t*30) for t in xrange(10)]
        random.seed(100)
        b = [rlc.get_flux(t*30) for t in xrange(10)]
        self.assertEqual(a, b)

    def test_02_output(self):
        rlc = RandomLightCurve()
        t1 = 10.0
        t2 = np.array(15)
        t3 = (10.0)
        #t4 = [10.0]
        a = rlc.get_flux(t1)
        b = np.array(rlc.get_flux(t2))
        c = rlc.get_flux(t3)
        #d = rlc.get_flux(t4)
        self.assertEqual(type(t1), type(a))
        self.assertEqual(type(t2), type(b))
        self.assertEqual(type(t3), type(c))
        #self.assertEqual(type(t4), type(d))

class AsciiLightCurveTestCase(TsigTestCase):

    def test_00_file_integrity(self):
        asciilc = ASCIILightCurve(name='tests/data/Kepler-6b.lc', cfgfile='tsig/lightcurve/example.cfg')
        asciilc.load_from_file(label='all')
        asciilc.write_to_file(outfile='tests/data/file_integrity_outfile.lc')
        filein = open('tests/data/Kepler-6b.lc')
        fileout = open('tests/data/file_integrity_outfile.lc')
        timein, timeout, magin, magout = [], [], [], []
        filein.next()
        fileout.next()
        for line in filein:
            timein.append(int(round(float(line.split()[0]), 7)*10**7))
            magin.append(int(round(float(line.split()[1]), 7)*10**7))
        for line in fileout:
            timeout.append(int(float(line.split()[0])*10**7))
            magout.append(int(float(line.split()[1])*10**7))
        self.assertEqual(timein, timeout)
        self.assertEqual(magin, magout)

    def test_01_num_cols(self):
        asciilc = ASCIILightCurve(name='tests/data/Kepler-6b_1col.lc', cfgfile='tsig/lightcurve/example.cfg')
        indexerror = False
        try:
            asciilc.load_from_file()
        except IndexError:
            indexerror = True
        self.assertTrue(indexerror)

    def test_02_get_flux_interpolation(self):
        asciilc = ASCIILightCurve(name='tests/data/Kepler-6b.lc', cfgfile='tsig/lightcurve/example.cfg')
        asciilc.load_from_file()
        time, flux = np.linspace(0, 1342, 61), [] #22 minute cadence
        for t in time:
            flux.append(asciilc.get_flux(t))
        self.assertEqual(len(time), len(flux))

    def test_03_input_exceed_lc_data(self):
        asciilc = ASCIILightCurve(name='tests/data/Kepler-6b.lc', cfgfile='tsig/lightcurve/example.cfg')
        asciilc.load_from_file()
        time, flux = np.linspace(0, 50000, 51), [] #test lightcurve only goes up to ~44,850 minutes
        for t in time:
            flux.append(asciilc.get_flux(t))
        self.assertEqual(len(time), len(flux))

class FITSLightCurveTestCase(TsigTestCase):

    def test_00_get_flux_interpolation(self):
        try:
            np.nanmedian([0])
        except AttributeError:
            print "skip test: numpy too old"
            return
        fitslc = FITSLightCurve(name='tests/data/fitstestfile.fits', cfgfile='tsig/lightcurve/example.cfg')
        fitslc.load_from_file()
        time, flux = np.linspace(0, 1342, 61), [] #22 minute cadence
        for t in time:
            flux.append(fitslc.get_flux(t))
        self.assertEqual(len(time), len(flux))

    def test_01_input_exceed_lc_data(self):
        try:
            np.nanmedian([0])
        except AttributeError:
            print "skip test: numpy too old"
            return
        fitslc = FITSLightCurve(name='tests/data/fitstestfile.fits', cfgfile='tsig/lightcurve/example.cfg')
        fitslc.load_from_file()
        time, flux = np.linspace(0, 50000, 51), [] #test lightcurve only goes up to ~44,850 minutes
        for t in time:
            flux.append(fitslc.get_flux(t))
        self.assertEqual(len(time), len(flux))

class ChelseaLightCurveTestCase(TsigTestCase):

    def test_00_choose_lc_from_file_extension(self):
        ascii_input = FilebasedLightCurve.create_lightcurve(infile='tests/data/Kepler-6b.lc', cfgfile='tsig/lightcurve/example.cfg')
        fits_input = FilebasedLightCurve.create_lightcurve(infile='tests/data/fitstestfile.fits', cfgfile='tsig/lightcurve/example.cfg')
        self.assertIsInstance(ascii_input, ASCIILightCurve)
        self.assertIsInstance(fits_input, FITSLightCurve)

class ReadColumnTestCase(TsigTestCase):

    def setUp(self):
        super(ReadColumnTestCase, self).setUp()
        self.asciilc = ASCIILightCurve(name='tests/data/Kepler-6b_erroneous.lc', cfgfile='tsig/lightcurve/example.cfg')
        self.data, self.cols, self.infile = [], 2, 'tests/data/Kepler-6b_erroneous.lc'

    def test_00_float_value_error_falsefill(self):
        with self.assertRaises(ValueError):
            self.asciilc.readcolumn(self.data, self.cols, self.infile, fill=False)

    def test_01_float_value_error_truefill(self):
        self.asciilc.readcolumn(self.data, self.cols, self.infile, fill=True)
        self.assertIs(self.data[3], np.nan)

    def test_02_float_index_error_falsefill(self):
        self.cols = 5
        with self.assertRaises(IndexError):
            self.asciilc.readcolumn(self.data, self.cols, self.infile, fill=False)

    def test_03_float_index_error_truefill(self):
        self.cols = 5
        self.asciilc.readcolumn(self.data, self.cols, self.infile, fill=True)
        self.assertIs(self.data[3], np.nan)

    def test_04_str_index_error_falsefill(self):
        self.cols = 5
        with self.assertRaises(IndexError):
            self.asciilc.readcolumn(self.data, self.cols, self.infile, datformat='str', fill=False)

    def test_05_str_index_error_truefill(self):
        self.cols = 5
        self.asciilc.readcolumn(self.data, self.cols, self.infile, datformat='str', fill=True)
        self.assertIs(self.data[3], None)

if __name__ == '__main__':
    numpy.random.seed(100) # ensure the same result each time
    test_support.run_unittest(LightCurveTestCase)
    test_support.run_unittest(RandomLightCurveTestCase)
    test_support.run_unittest(AsciiLightCurveTestCase)
    test_support.run_unittest(FITSLightCurveTestCase)
    test_support.run_unittest(ChelseaLightCurveTestCase)
    test_support.run_unittest(ReadColumnTestCase)
