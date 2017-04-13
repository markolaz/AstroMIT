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
import optparse

from lctools.util.util import getext 
from lctools.util.configurable import ConfigurableObject

# create logger
logger = logging.getLogger(__name__)



class VtDft(ConfigurableObject):
    LENGTH = 3 
    config_keys = ['DFT']
    def __init__(self, nbeam=5, f1=10., npeaks=3, nclip=5, nclipiter=1, outdir='', analfile='.dftanal'):       
        super(VtDft, self).__init__()
        self.nbeam = nbeam	
        self.f1 = f1		
        self.npeaks = npeaks 
        self.nclip = nclip 
        self.nclipiter = nclipiter 
        self.outdir = outdir
        self.analfile = analfile
        logger.debug("nbeam = %d, f1=%f, npeaks=%d, nclip=%d, nclipiter=%d, outdir=%s, analfile=%s", self.nbeam, self.f1, self.npeaks, self.nclip, self.nclipiter, self.outdir, self.analfile)
    	return
    def __call__(self, lcfile, replace=False):
        if self.outdir == '':
            outdir = os.path.dirname(lcfile.name) + '/'
        else:
            outdir = self.outdir
        dftanalf = outdir + self.analfile

        dftanalfile = getext(lcfile.name, dftanalf, '')
        dftpath = os.path.dirname(dftanalfile)
        logger.debug("dftanalfile=%s, dftpath=%s", dftanalfile, dftpath)
        if dftpath == "":
            dftpath = "./"
        modelpath = dftpath
        phasepath = dftpath
        if os.path.exists(dftanalfile) and (not replace):
            logger.warning("File %s exists, do not allow to replace, set -r to force overwrite", dftanal)
            return
        elif (not os.path.exists(lcfile.name)):
            logger.warning("Input file %s do not exists, do nothing", lcfile.name)
            return

        else:
            cmdline = "vartools -i %s -header -readformat 1 %d %d 3 -dftclean %d maxfreq %f finddirtypeaks %d clip %d %d clean 0.5 5.0 outcspec %s findcleanpeaks %d clip %d %d verboseout" % (lcfile.name, lcfile.cols['jd'], lcfile.cols['ltflc'], self.nbeam, self.f1, self.npeaks, self.nclip, self.nclipiter, dftpath, self.npeaks, self.nclip, self.nclipiter)
            logger.info("excuting cmdline: %s", cmdline)
            status, output = commands.getstatusoutput(cmdline)
            # print output
            header, result = output.split('\n')
            newheader = " ".join(header.split()[1:VtDft.LENGTH+1])
            newlines = "" 
            for i in xrange(self.npeaks):
                newlines += "%d " % (i+1) 
                newlines += " ".join(result.split()[1+i*VtDft.LENGTH:(1+VtDft.LENGTH*(i+1))])+"\n"
            fout = open(dftanalfile, mode="w")
            fout.write("#")
            fout.write(newheader)
            fout.write("\n")
            fout.write(newlines.rstrip())
            fout.close()
        return 


