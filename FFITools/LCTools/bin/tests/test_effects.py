#!/usr/bin/env python
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
Test tsig.effects API in a general way.
"""

import os
import sys

from base import TsigTestCase, test_support

from tsig.effects.base import Effect, TestEffect

class EffectsApiTestCase(TsigTestCase):
    """Unit tests for all classes in the tsig.effects module.
    
    Like the camera test we test importing known effects, but also
    any unit tests of the base class should go here.
    """
    def test_01_process(self):
        effect = TestEffect()
        effect.feed_in(self.assertDataFile('example.fits'))
        effect.apply_now()
        effect.feed_out(True, 'test')
        self.assertDataFile('example__test__00__testeffect.fits')

    def test_02_errors(self):
        """Test basic errors"""
        effect = Effect()
        with self.assertRaises(ValueError) as context:
            effect.apply_now()
        effect.feed_in(self.assertDataFile('example.fits'))
        with self.assertRaises(NotImplementedError) as context:
            effect.apply_now()
        with self.assertRaises(IOError) as context:
            effect.feed_out(True)



if __name__ == '__main__':
    test_support.run_unittest(EffectsApiTestCase)

