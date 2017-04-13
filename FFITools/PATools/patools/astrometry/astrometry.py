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




import numpy as np
import scipy as sp
import os
import logging
from patools.util.configurable import ConfigurableObject
from starlist import Starlist
# create logger
logger = logging.getLogger(__name__)


class Astrometry(object):

    def __init__(self):
        return

    def __call__(self, infile, outfile=''):
        return

# FIXME: configurable object problem, see fistar.py

class Grmatch(ConfigurableObject):
    config_keys = ['Grmatch', 'Catalog']
    def __init__(self, order=2, maxdistance=1, unitarity=0.01, catfile=None):
        super(Grmatch, self).__init__()
        self.order = int(order)
        self.maxdistance = int(maxdistance)
        self.unitarity = float(unitarity)
        self.catfile = catfile
        logger.debug("order=%d, maxdistance=%d, unitarity=%f, catfile=%s", self.order, self.maxdistance, self.unitarity, self.catfile)
        if self.catfile is None or (not os.path.exists(self.catfile)):
            raise IOError, 'a catalog file is required for ' \
                           'the astrometry step, file %s does ' \
                           'not exists' % self.catfile
        
        self.catalog = Starlist(self.catfile)
    def get_transfile(self, infile): 
        transfile = os.path.splitext(infile)[0] + '.trans'
        return transfile

    def proj_sky_to_xy(self, xylist, inputcat='', transfile='', xlim=2048, ylim=2048):
        if transfile == '':
            transfile = self.get_transfile(xylist)
        if inputcat == '':
            inputcat = self.catalog
        cmd_line = "grtrans --comment --input %s --col-xy %d,%d " \
                 "--input-transformation %s --col-out %d,%d --output - | awk '$%d>0 && $%d<%d && $%d>0 && $%d<%d {print}' >%s" \
                 % (inputcat.name, inputcat.colx, inputcat.coly,
                    transfile, inputcat.colx, inputcat.coly, inputcat.colx, inputcat.colx, xlim, inputcat.coly, inputcat.coly, ylim, xylist)
        logging.debug("Excute projection: %s", cmd_line)
        os.system(cmd_line)
        return


    def __call__(self, starlist, outfile='', matchfile=''):
        if outfile == '':
            outfile = self.get_transfile(starlist.name) 
        cmdline = 'grmatch --match-points -r ' + self.catalog.name \
                  + ' --col-ref %d,%d --col-ref-ordering -%d -i ' \
                    % (self.catalog.colx, self.catalog.coly, self.catalog.colmag) \
                  + starlist.name \
                  + ' --col-inp %d,%d --col-inp-ordering +%d ' \
                    '--weight reference,column=%d ' \
                    '--triangulation maxinp=1000,maxref=1000,auto,unitarity=%.2f ' \
                    '--order %d --max-distance %2.1f --comment --output-transformation %s ' \
                    % (starlist.colx, starlist.coly, starlist.colmag,
                       self.catalog.colmag, self.unitarity, self.order,
                       self.maxdistance, outfile)
        if not matchfile == '':
            cmdline += '--output-matched %s ' % matchfile
        logging.debug("Excute grmatch: %s", cmdline)
        os.system(cmdline)
        return


