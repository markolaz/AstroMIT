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


import logging
import os
from patools.star import Fistar
from patools.astrometry import Grmatch, Starlist, Catalog
from patools.phot import Fiphot

# create logger
logger = logging.getLogger(__name__)


class Frame(object):
    def __init__(self, cfg):
        self.cfg = cfg 
        self.base = ''
        self.indir = self.cfg.get('Setup').get('indir')
        self.fistar = Fistar(self.cfg)
        self.fiphot = Fiphot(self.cfg)
        self.catalog = Catalog(self.cfg)
        self.catbase = self.catalog.name
        self.catalog.name = os.path.splitext(self.catbase)[0]+'_bright.txt'
        if not os.path.exists(self.catalog.name):
            logger.info("Query small catalog for bright stars")
            self.catalog.query(maglim=11)
        self.catalog.project()
        self.catalog.name = os.path.splitext(self.catbase)[0]+'_full.txt'
        if not os.path.exists(self.catalog.name):
            logger.info("Query large catalog for all stars")
            self.catalog.query(maglim=14)
            self.catalog.project()

        self.solver = Grmatch(self.cfg)

        self.catproj = None
        
        astrom_arg = None
        if not astrom_arg is None: 
            self.catname = eval(p.get('Catalog', 'catfile'))
            
            self.ra0, self.dec0, self.ra, self.dec, self.width = self.getheader()
                
        return
    def get_filenames(self):
        if not self.base =='':
            self.fitsfile = self.indir + self.base + '.fits'
            self.xylist = self.indir + self.base + '.xyls'
            self.transfile = self.indir + self.base + '.trans'
            self.fistarfile = self.indir + self.base + '.fistar'

    def getheader(self):
        # need to somehow figure out the field from fitsheader
        p = ConfigParser.RawConfigParser()
        p.read(self.cfgfile)
        # hack for the simple example
        try:
            ra0 = eval(p.get('Catalog', 'ra0'))
        except:
            ra0 = 270.0
        try:
            dec0 = eval(p.get('Catalog', 'dec0'))
        except:
            dec0 = 66.56
        try: 
            ra = eval(p.get('Catalog', 'ra'))
        except:
            ra = 249.36
        try:
            dec = eval(p.get('Catalog', 'dec'))
        except:
            dec = 71.86
        try:
            width = eval(p.get('Catalog', 'width'))
        except:
            width = 9.

        return [ra0, dec0, ra, dec, width]
    
    

    def extract(self, dry=False):
        logger.info("Begining source extraction")
        self.fistar(self.fitsfile)
        starlist = Starlist(self.fistarfile, colx=2, coly=3, colmag=9) 
        logger.info("Begining solve for astrometry")
        self.solver(starlist)
        self.solver.proj_sky_to_xy(self.xylist)
        catproj = Starlist(self.xylist)
        logger.info("Begining photometry")
        self.fiphot(self.fitsfile, catproj)
        logger.info("Finish extracting photometry from frame")

        return


