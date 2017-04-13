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


import os
import time
import tempfile
import commands
import logging
import sys

from lctools.util.util import getext 
from lctools.util.configurable import ConfigurableObject

# create logger
logger = logging.getLogger(__name__)

class VtBls(ConfigurableObject):
    LENGTH = 20
    config_keys = ["BLS"]
    """
    A python class that configure and run the vartools BLS algorithm.
    """
    def __init__(self, f0=0.0625, f1=2.0, fn=100000, qmin=0.008, qmax=0.08, nbin=200, peaknum=3, outdir='', analfile='.blsanal'):
        super(VtBls, self).__init__()
        self.f0 = float(f0)
        self.f1 = float(f1)
        self.fn = float(fn)
        self.qmin = float(qmin)
        self.qmax = float(qmax)
        self.nbin = float(nbin)
        self.peaknum = int(peaknum)
        self.outdir = outdir
        self.analfile = analfile
        logger.debug("readfrom configure file, f0=%f, f1=%f, fn=%d, qmin=%f, qmax=%f, nbin=%d, peaknum=%d, outdir=%s, analfile=%s", self.f0, self.f1, self.fn, self.qmin, self.qmax, self.nbin, self.peaknum, self.outdir, self.analfile)
        return
    def get_blsanalfile(self, infile):
        blsanalf = self.outdir+self.analfile
        blsanalfile = getext(infile, blsanalf, '')
        return blsanalfile
    def __call__(self, lcfile, replace=False):
        if self.outdir == '':
            outdir = os.path.dirname(lcfile.name) + '/'
        else:
            outdir = self.outdir
        blsanalf = outdir + self.analfile
        blsanalfile = getext(lcfile.name, blsanalf, '')
        blspath = os.path.dirname(blsanalfile)
        logger.debug("outdir=%s, blsanalf=%s, blsanalfile=%s, blspath=%s", outdir, blsanalf, blsanalfile, blspath)
        if blspath == "":
            blspath = "./"
        modelpath = blspath
        phasepath = blspath
        if os.path.exists(blsanalfile) and (not replace):
            logger.warning("File %s exists, do not allow to replace, set -r to force overwrite", blsanalfile)
            return
        elif (not os.path.exists(lcfile.name)):
            logger.warning("Input file %s do not exists, do nothing", lcfile.name)
            return

        else:
            cmdline = "vartools -i %s -header " \
                      "-readformat 1 %d %d 3 " \
                      "-BLS q %f %f %f %f %d %d %d %d 1 %s 1 %s 0 " \
                      "fittrap nobinnedrms ophcurve %s -0.1 1.1 0.001" \
                      % (lcfile.name, lcfile.cols['jd'], lcfile.cols['ltflc'],
                         self.qmin, self.qmax, 1./self.f1, 1./self.f0, self.fn,
                         self.nbin, 0, self.peaknum, blspath, modelpath, phasepath)
            logger.info("excuting command line %s", cmdline)
            status, output = commands.getstatusoutput(cmdline)
            print output
            header, result = output.split('\n')
            newheader = " ".join(header.split()[1:VtBls.LENGTH+1])
            newlines = ""
            for i in xrange(self.peaknum):
                newlines += "%d " % (i + 1)
                newlines += " ".join(result.split()[1 + i * VtBls.LENGTH:(1 + VtBls.LENGTH * (i + 1))]) + "\n"
            fout = open(blsanalfile, mode="w")
            fout.write("# BLS_No ")
            fout.write(newheader)
            fout.write("\n")
            fout.write(newlines.rstrip())
            fout.close()
            return
        

