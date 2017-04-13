from setuptools import setup

import os

if os.path.exists('MANIFEST'): os.remove('MANIFEST')

with open('README') as fhl:
    description = fhl.read()

setup(
        name = __pkgname__,
        version = __version__,
        description = 'FFITools', 
        long_description = description,
        author = 'Tess Science',
        url = ''
        author_email = '' 
        test_suite = 'tests',
        platforms = 'linux',
        license = 'GPLv3',
        py_modules = ['PATools', 'LCTools', 'DVProducts',],
        classifiers = [
            ],
        )

