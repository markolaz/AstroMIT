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
The effects base classes provides all the tools required to write
a photon or electron based effect. These should be applied to either
a Camera() object or a CCD() object depending on what you need.
"""

import re
import os
import sys
import time
import inspect

import numpy as np
from astropy.io import fits

class classproperty(object):
    def __init__(self, f):
        self.f = f

    def __get__(self, obj, owner):
        return self.f(owner)


class Effect(object):
    """
    The base class for any effect which effects fits data directly.
    """
    priority = None

    def __init__(self):
        """
        Initiaize the effect, your child class should implement it's own
        init which calls this one as 'super' before saving it's own params.
        """
        self.filename = None
        self.fits = None

    def apply(self, hdulist):
        """
        The apply function is where the work of the effect happens. Your child
        class should overwrite this and the effect will be run.

          hdulist - The fits data, passed through for convenience.

        """
        raise NotImplementedError("An effect must have an apply() function.")

    def feed_in(self, fits_in):
        """
        Used by Camera and CCD to chain together effects, do not override
        unless you know what you're doing.

          fits - A fits filename or a HDUList fits file.

        """
        if isinstance(fits_in, basestring):
            if os.path.isfile(fits_in):
                self.filename = fits_in
                self.hdulist = fits.open(fits_in, memmap=True)
            else:
                raise IOError("Can't load fits file: %s" % fits_in)
        elif isinstance(fits_in, fits.HDUList):
            self.hdulist = fits_in
        else:
            raise ValueError("Not sure what type of data this "
                    "is: %s" % type(fits_in).__name__)

    def feed_out(self, to_file=False, run=None, overwrite=True):
        """
        Used by Camera and CCD to chain together feeds, do not override
        unless you know what you're doing.

          to_file
            - If true, the fits data is written to a filename before being
              returned as the new filename for the next effect.
            - If false, the HDUList is returned to the chain.

        """
        if not hasattr(self, 'hdu_output'):
            raise IOError("Effect has not been applied yet so can not output.")

        if to_file:
            if self.filename is None:
                raise IOError("You must feed a filename into the effect to get"
                        " a fits file out.")
            template = '%(filename)s__%(run)s__%(index)02d__%(name)s.%(ext)s'
            (fn, ext) = self.filename.rsplit('.', 1)
            if '__' in fn:
                (fn, run, index, _) = fn.split('__')
                index = int(index) + 1
            else:
                run = run or str(int(time.time()))
                index = 0
            new_filename = template % {
              'filename': fn,
              'run': run,
              'index': index,
              'name': self.name.lower(),
              'ext': ext,
            }
            self.hdu_output.writeto(new_filename, overwrite=overwrite)
            return new_filename
        return self.hdu_output

    def apply_now(self):
        """
        Internally run the apply function and save the returned hdulist
        data to hdu_output which is then used by out() and saved correctly.
        """
        if not hasattr(self, 'hdulist'):
            raise ValueError("feed_in must be called before apply.")

        self.hdu_output = self.apply(self.hdulist)

        self.hdu_output[0].header['comment'] = 'Applied Effect: %s' % self.name

    @classproperty
    def name(cls):
        """By default returns the class name as the name of the effect."""
        return cls.__name__

    @classproperty
    def title(cls):
        return re.sub("([a-z])([A-Z])","\g<1> \g<2>", cls.name)

    @classproperty
    def args(cls):
        """Return a list of arguments and their default values"""
        docs = (cls.__init__.__doc__ or "").split("\n")
        docs = (line.strip().split(" - ", 1) for line in docs)
        docs = dict((l[0].strip(), l[1].strip()) for l in docs if len(l) == 2)

        spec = inspect.getargspec(cls.__init__)
        args = spec.args[1:]
        defs = spec.defaults or ()
        defaults = ((None,) * (len(args) - len(defs))) + defs

        for x, arg in enumerate(args):
            default = defaults[x]
            if default is not None:
                yield (arg, default, type(default), docs.get(arg, None))
            else:
                yield (arg, None, str, docs.get(arg, None))

    @classproperty
    def all(cls):
	subclasses = set()
	work = [cls]
	while work:
	    parent = work.pop()
	    for child in parent.__subclasses__():
		if child not in subclasses:
		    subclasses.add(child)
		    work.append(child)
	return list(subclasses)


class ImageEffect(Effect):
    def apply_to_image(self, image):
        """
        The main entry point to applying an effect to an image.

          image - Is a np.array converted to float64 from int

        """
        raise NotImplementedError("apply_to_image() has not been implemented")

    def apply(self, hdulist):
        """
        Performs the conversion to and from images formats as needed
        """
        image = hdulist[0].data.astype(float)
        image = self.apply_to_image(image)
        hdulist[0].data = image.astype(np.int32)
        for (name, value) in self.get_headers():
            hdulist[0].header.append((name, value[0], value[1]))
        return hdulist

    def get_headers(self):
        return []


class TestEffect(Effect):
    def apply(self, data):
        return data

