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
# but ITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License 
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 



import numpy as np
import scipy as sp
import os
import sys
import commands
import tempfile


from patools.frame import Frame 
from patools.util.configurable import ConfigurableObject

# create logger
logger = logging.getLogger(__name__)



class Field(Configurableobject):
    def __init__(self, indir='', nccds=4, ncameras=4, nexps=1248, basename='tsig', ra0=None, dec0=None):
        print "reading paramters from %s" % cfgfile
        self.frames = []
        self.basedir = indir 
        self.nccds = nccds 
        self.nexps = nexps 
        self.ncameras = ncameras 
        self.basename = basename 

        # init field center ra,dec
        self.ra0 = ra0 
        self.dec0 = dec0
        if not ra0 or not dec0
            self.get_header()

        # init ccd center ra,dec, ccdwidth
        self.ccdra = np.zeros(self.nccds)
        self.ccddec = np.zeros(self.nccds)
        self.width = 0
        self.get_ccd_center()

    def get_header(self):
        # need to be more flexible
        frame0 = self.basepath+"ccd1/"+self.basename+"_ccd1_%.6d.fits" % 0
        cmd_line = "fiheader --read CRVAL1,CRVAL2,SCALE,NAXIS1 %s" % frame0
        status, output = commands.getstatusoutput(cmd_line)
        print output
        self.ra0 = eval(output.split('\n')[0].split()[1])
        self.dec0 = eval(output.split('\n')[1].split()[1])
        scale = eval(output.split('\n')[2].split()[1])
        nx = eval(output.split('\n')[3].split()[1])
        self.width = scale*nx/3600.*np.sqrt(2.)
        return 
    
    def get_ccd_center(self):
        # currently we resort to astormetry.net for a first solution, might
        # need to change later on we may also save the wcs solution somewhere
        # so we don't repeat this procedure every time
        def get_ra_dec_anet(frame):
            tempwcs = os.path.splitext(frame)[0]+".wcs"
            if not os.path.exists(tempwcs):
                cmd_line = "solve-field --overwrite %s --ra %f --dec %f " \
                           "--radius %d -p -q 0.01 --tweak 5 --wcs %s" \
                           % (frame, self.ra0, self.dec0, 10, tempwcs)
                os.system(cmd_line)
            print tempwcs 
            status, result = commands.getstatusoutput("wcsinfo %s | grep 'ra_center '" % tempwcs)
            ra = eval(result.split()[1])
            status, result = commands.getstatusoutput("wcsinfo %s | grep 'dec_center '" % tempwcs)
            dec = eval(result.split()[1])
            return [ra, dec]

        for i in xrange(self.nccds):
            frame0 = self.basepath+"ccd%d/" % (i+1)+self.basename+"_ccd%d_%.6d.fits" % (i+1, 0)
            self.ccdra[i], self.ccddec[i] = get_ra_dec_anet(frame0)
        return

    def init_frames(self):
        for c in xrange(self.nccds):
            inpath = self.basepath+"ccd%d/" % (c+1)
            catname = self.basepath+"ccd%d/catalog.txt" % (c+1)
            astrom_arg = [catname, self.ra0, self.dec0, self.ccdra[c], self.ccddec[c], self.width]
            for i in xrange(self.ncadences):
                framebase = self.basename+"_ccd%d_%.6d" % (c+1, i)
                self.frames.append(Frame(framebase, cfgfile=self.cfgfile, inpath=inpath, astrom_arg=astrom_arg))
        return

    def extract(self, dry=False):
        for frame in self.frames: 
            frame.extract(dry=dry)
            print "Done with Frame %s" % frame.base


