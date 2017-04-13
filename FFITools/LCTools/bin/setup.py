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

from setuptools import setup
from tsig import __version__, __pkgname__
import os

# remove MANIFEST. distutils doesn't properly update it when the
# contents of directories change.
if os.path.exists('MANIFEST'): os.remove('MANIFEST')

# Grab description for Pypi
with open('README.md') as fhl:
    description = fhl.read()

bin_dir = os.path.join(os.path.dirname(__file__), 'bin')

setup(
    name             = __pkgname__,
    version          = __version__,
    description      = 'TSIG',
    long_description = description,
    author           = 'TESS Science',
    url              = 'https://tessgit.mit.edu/tess/tsig',
    author_email     = 'tess-info@mit.edu',
    test_suite       = 'tests',
    platforms        = 'any',
    license          = 'GPLv3',
    scripts=[os.path.join('bin', a) for a in os.listdir(bin_dir)],
    requires=[
        'configobj(>=5.0.6)', # configuration files
        'Pillow(>=2.6.1)', # used for test patterns
        'astropy(>=1.3)',
        'matplotlib(>=2.0.0)',
        'numpy(>=1.12.0)',
        'scipy(>=0.14.0)',
        'psycopg2(>=2.5.4)', # TIC queries
    ],
    packages=[
        'tsig', 'tsig.camera', 'tsig.catalog', 'tsig.effects',
        'tsig.lightcurve', 'tsig.spacecraft', 'tsig.util',
    ],
    classifiers=[
      'Development Status :: 1 - Planning',
      'Intended Audience :: Developers',
      'Intended Audience :: Information Technology',
      'Intended Audience :: Science/Research',
      'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
      'Operating System :: POSIX',
      'Operating System :: POSIX :: Linux',
      'Programming Language :: Python',
      'Programming Language :: Python :: 2.7',
      'Programming Language :: Python :: 3.3',
    ],
    include_package_data=True,
 )
