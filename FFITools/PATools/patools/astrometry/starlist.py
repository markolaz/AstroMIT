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




from astroquery.vizier import Vizier
import astropy.coordinates as coord
import astropy.units as u
import tempfile
import os
import numpy as np
from patools.util.dataio import readcolumn
from patools.util.configurable import ConfigurableObject
import logging
logger = logging.getLogger(__name__)

class Starlist(object):

    def __init__(self, infile, colid=1, colra=2, coldec=3, colmag=4,
                 colx=5, coly=6):
        self.name = infile
        self.colid = colid
        self.colra = colra
        self.coldec = coldec
        self.colmag = colmag
        self.colx = colx
        self.coly = coly
        return


class Catalog(ConfigurableObject):
    config_keys = ['Catalog']
    def __init__(self, catfile='', ra=249.36, dec=71.86, width=13, colid=1, colra=2, coldec=3,
                 colmag=4, colx=5, coly=6, ra0=270.0, dec0=66.56, method=None):
        self.name = catfile
        self.ra = float(ra)
        self.dec = float(dec)
        self.ra0 = float(ra0)
        self.dec0 = float(dec0)
        if self.ra0 is None:
            self.ra0 = self.ra
        if self.dec0 is None:
            self.dec0 = self.dec
        self.width = float(width)
        self.colid = colid
        self.colra = colra
        self.coldec = coldec
        self.colmag = colmag
        self.colx = colx
        self.coly = coly
        self.method = method

        logging.debug("Config Catalog: name=%s, ra=%f, dec=%f, ra0=%f, dec0=%f, width=%f, method=%s", self.name, self.ra, self.dec, self.ra0, self.dec0, self.width, self.method)
    def query(self, maglim=11.):
        if self.method is None:
            # need to replace this with lfov-gaia or lfov-usno
            v = Vizier(columns=['UCAC4', '_RAJ2000', 'e_RAJ2000', '_DEJ2000',
                                'e_DEJ2000', 'Jmag', 'Kmag', 'Vmag', 'pmRA', 'pmDE',
                                'e_pmRA', 'e_pmDE', 'rmag', 'imag'],
                       row_limit="unlimited", column_filters={"Vmag": "<%f" % maglim})
            # print self.ra,self.dec
            result = v.query_region(coord.SkyCoord(ra=self.ra*u.degree,
                                                   dec=self.dec*u.degree,
                                                   frame='icrs'),
                                    width=self.width*u.degree, catalog=["UCAC4"])

            ids = result[0]['UCAC4']
            ras = result[0]['_RAJ2000'][:]
            decs = result[0]['_DEJ2000'][:]
            e_ra = result[0]['e_RAJ2000'][:]
            e_de = result[0]['e_DEJ2000'][:]
            jmag = result[0]['Jmag'][:]
            kmag = result[0]['Kmag'][:]
            vmag = result[0]['Vmag'][:]
            fout = open(self.name, mode='w')
            for i in xrange(len(ids)):
                # ID[i]=''.join(ID[i].split('-'))
                if jmag[i] == "--":
                    jmag[i] = np.nan
                if kmag[i] == "--":
                    kmag[i] = np.nan
                if vmag[i] == "--":
                    vmag[i] = np.nan
                fout.write("%s %12.7f %12.7f %5.3f %d %d %5.3f %5.3f %5.3f\n" %
                           (ids[i], ras[i], decs[i], vmag[i], e_ra[i], e_de[i],
                            jmag[i], kmag[i], vmag[i]))
            fout.close()
        elif self.method == 'lfov-usno':
            cmd_line = "lfov-usno.sh %25.12f %25.12f %12.7f > %s" % \
                     (self.ra, self.dec, self.width, self.name)
            os.system(cmd_line)
            print cmd_line
            # raise AttributeError, "not implemented yet"
        else:
            raise AttributeError, "not implemented yet for catalog " \
                                  "query method %s" % self.method
        return

    

    def project(self):
        tmpfile = tempfile.NamedTemporaryFile(delete=False)
        cmd_line = 'grtrans --input %s --wcs tan,degrees,ra=%.12f,dec=%.12f ' \
                   '--col-radec 2,3 --col-out 5,6 --output %s' % \
                   (self.name, self.ra0, self.dec0, tmpfile.name)
        logger.debug("project catalog: %s", cmd_line)
        os.system(cmd_line)
        ids = []
        readcolumn(ids, 1, tmpfile.name, datformat='str')
        ids = np.array(ids)
        data = np.loadtxt(tmpfile.name, usecols=(1, 2, 3, 4, 5))
        x0 = np.cos(self.dec * np.pi/180.) * np.cos(self.ra * np.pi/180.)
        y0 = np.cos(self.dec * np.pi/180.) * np.sin(self.ra * np.pi/180.)
        z0 = np.sin(self.dec * np.pi/180.)
        zmin = np.cos(self.width * np.pi/180.)
        # x=data[:,3]
        # y=data[:,4]
        p_ra = data[:, 0]
        p_dec = data[:, 1]
        cx = np.cos(p_dec * np.pi/180.) * np.cos(p_ra * np.pi/180.)
        cy = np.cos(p_dec * np.pi/180.) * np.sin(p_ra * np.pi/180.)
        cz = np.sin(p_dec * np.pi/180.)
        w = cx * x0 + cy * y0 + cz * z0
        index = zmin <= w
        newdata = data[index, :]
        newids = ids[index]
        fout = open(self.name, mode='w')
        fout.write("#ID RA DEC Mag Chi Eta ...\n")
        fout.write("#[1] [2] [3] [4] [5] [6] ...\n")
        for i in xrange(newdata.shape[0]):
            fout.write("%s " % newids[i])
            for j in xrange(newdata.shape[1]):
                fout.write("%12.7f " % newdata[i, j])
            fout.write("\n")
        fout.close()
        os.unlink(tmpfile.name)
        return
