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


import math
import numpy as np
import scipy as sp
import copy
from scipy import linalg
from scipy.optimize import leastsq
from scipy import interpolate
from scipy import signal
from scipy import stats
import ConfigParser
import logging

import lctools 
from lctools.util.configurable import ConfigurableObject

# create logger
logger = logging.getLogger(__name__)

class DetrendFunc(ConfigurableObject):
    config_keys = ['LC']
    def __init__(self):
        super(DetrendFunc, self).__init__()
        self.required_cols = ['jd', 'rlc']
        return
    
    def __getstate__(self):
        d = self.__dict__.copy()
        if 'logger' in d.keys():
            d['logger'] = d['logger'].name
        return d

    def __setstate__(self, d):
        if 'logger' in d.keys():
            d['logger'] = logging.getLogger(d['logger'])
        self.__dict__.update(d)

   
    @staticmethod
    def create(cfg, method='cos'): 
        for c in [PolyDetrend, COSDetrend, EPDDetrend, FocusDetrend]:
            if c.istype(method):
                return c(cfg)

        raise NotImplementedError, 'detrending method %s is not recognized' % method 
    def check_data(self, lightcurve):
        for key in self.required_cols:
            lightcurve.check_data(key)
        return 

    def detrend(self, lightcurve):
        self.check_data(lightcurve)
        ltf = self._detrend(lightcurve.data)
        lightcurve.data['ltflc'] = ltf
        return  
    
    def _detrend(data):
        raise NotImplementedError
        

class PolyDetrend(DetrendFunc):
    # Detrend with polynomial fit
    METHOD = 'poly'
    def __init__(self, wn=13, norder=5):
        super(PolyDetrend, self).__init__()
        self.wn = wn
        self.norder = norder
        
        self.required_keys = ['wn', 'norder']

    @staticmethod
    def istype(method):
        return method == PolyDetrend.METHOD


    def _detrend(self, data, noplot=True, intran=[]):
        otime = data['jd']
        if not intran:
            intran = otime < 0
        oflux = data['rlc']
        order = self.norder
        length = len(oflux)
        ctime = otime[-intran]
        cflux = sp.signal.medfilt(oflux[-intran] - np.mean(oflux[-intran]), self.wn)
        fitflux, c = lspolyordern(ctime, cflux, order)
        rflux = np.zeros(length)
        for i in range(order + 1):
            rflux += c[i]*otime**(order-i)
        
        if not noplot:

            import matplotlib
            import matplotlib.pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(otime, oflux, '.')
            plt.show()

            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(otime, rflux, '.', ctime, cflux, 'x')
            plt.show()
        dflux = oflux - rflux
        dflux -= np.mean(dflux) - np.mean(oflux)
        
        return dflux


class COSDetrend(DetrendFunc):
    # Detrend with high pass filter
    METHOD = 'cos'
    def __init__(self, tmin=1.0, wn=13, ncomp=30):
        super(COSDetrend, self).__init__()
        self.required_keys = ['tmin', 'wn', 'ncomp']
        self.tmin = float(tmin)
        self.wn = int(wn)
        self.ncomp = int(ncomp)
        logger.info("USE a high pass filter to detrend, it is configured with %s", self.__str__)
    
    def __str__(self):
        parastr=''
        for key in self.required_keys:
            parastr+="%s=%s, " % (key, getattr(self,key))
        return parastr

    @staticmethod
    def istype(method):
        return method == COSDetrend.METHOD

    def _detrend(self, data, noplot=True, intran=[]):
        otime = data['jd'] 
        oflux = data['rlc']
        length = len(oflux)
        index = (~np.isnan(oflux))
        ctime = otime[index]
        cflux = oflux[index] - np.nanmean(oflux[index])
        cflux = sp.signal.medfilt(cflux, self.wn)

        if not noplot:
            import matplotlib
            from matplotlib import pyplot as plt
            plt.plot(otime, oflux-np.nanmean(oflux[index]), '+')
            plt.plot(ctime, cflux, 'r.')
            plt.xlim(7740,7750)
            plt.show()
        k = (cflux[-1] - cflux[0]) / (ctime[-1] - ctime[0])
        b = cflux[0] - k * ctime[0]
        cflux -= (k * ctime + b)
        e0 = min(ctime)
        timespan = max(ctime) - min(ctime)
        logger.debug("E0=%f, len(cflux)=%d, len(ctime)=%d", e0, len(cflux),len(ctime))
        amatrix = matrixgen(ctime-e0, self.ncomp, timespan)
        # print A.shape,n
        c, resid, rank, sigma = linalg.lstsq(amatrix, cflux)
        n = self.ncomp
        # print resid
        e = np.rollaxis(np.sin(np.array(np.r_['c', 0: n]
                                        * (otime[np.arange(length)]-e0))
                               * math.pi/timespan), 1, 0)
        rflux = np.dot(e, c)
        eflux = k * otime + b
        # print rflux.shape,eflux.shape
        rflux += eflux
        if not noplot:

            import matplotlib
            import matplotlib.pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(otime, rflux, '.', otime, eflux, '+', ctime, cflux, 'x')
            plt.show()
        dflux = oflux - rflux
        dflux -= np.nanmean(dflux)
        dflux += np.nanmean(oflux)
        return dflux


class EPDDetrend(DetrendFunc):
    # Detrend with external parameters
    METHOD = 'epd'
    def __init__(self):
        super(EPDDetrend, self).__init__()
        self.required_cols = ['jd', 'rlc', 'x', 'y', 'bg', 'bgerr']
        return
    
    @staticmethod
    def istype(method):
        return method == EPDDetrend.METHOD


    def _detrend(self, data):
        centx = data['x']
        centy = data['y']
        bg = data['bg']
        bgerr = data['bgerr']
        mag = data['ltf']

        def e(c, x, y, bg, bgdev, mag):
            return c[0] + c[1] * np.sin(2 * np.pi * x) \
                   + c[2] * np.cos(2 * np.pi * x) \
                   + c[3] * np.sin(2 * np.pi * y) \
                   + c[4] * np.cos(2 * np.pi * y) \
                   + c[5] * np.sin(4 * np.pi * x) \
                   + c[6] * np.cos(4 * np.pi * x) \
                   + c[7] * np.sin(4 * np.pi * y) \
                   + c[8] * np.cos(4 * np.pi * y) \
                   + c[9] * bg \
                   + c[10] * bgdev \
                   - mag
        # Kepmag = eval(result)
        c0 = np.zeros(11) + 1.
        c1, success = leastsq(e, c0, args=(centx, centy, bg, bgerr, magarr[a, :]))
        if (c1 == 1).all():
            dmag = mag * 1.0
            logger.warning("epd failed") 
        else:
            dmag = -e(c1, centx, centy, bg, bgerr, mag)
        dmag = dmag - sp.stats.nanmedian(dmag) + sp.stats.nanmedian(mag)
        return dmag


class FocusDetrend(DetrendFunc):
    # Detrend with Focus series only
    METHOD = 'focus'
    def __init__(self, wn=13, norder=3):
        super(FocusDetrend, self).__init__()
        self.required_cols = ['cadence', 'rlc']
        self.required_keys = ['wn', 'norder']
        self.norder = norder
        self.wn = wn
        if not cfgfile == '':
            self.config(cfgfile)
        return
    
    @staticmethod
    def istype(method):
        return method == FocusDetrend.METHOD



    def _detrend(self, data):
        xf = ((data['cadence'] / 48. / 13.7 + 0.5) % 1) - 0.5
        mag = data['rlc']
        focus = 0.0 + 10.0 * np.exp(-0.5 * (xf / 0.1)**2.)
        dmag = np.zeros(len(mag)) + mag
        noflux = sp.signal.medfilt(mag, self.wn)
        n = self.norder
        seg1 = data['cadence'] < 624
        x = data['cadence'][seg1]
        xi = sp.zeros(len(x)) + 1
        xii = sp.zeros(len(x)) + focus[seg1]
        amatrix = (x**n)[:, np.newaxis]
        for i in xrange(1, n):
            amatrix = np.hstack((amatrix, (x**(n-i))[:, np.newaxis]))
        amatrix = np.c_[amatrix, xi[:, np.newaxis]]
        amatrix = np.c_[amatrix, xii[:, np.newaxis]]
        c, resid, rank, sigma = sp.linalg.lstsq(amatrix, noflux[seg1])
        z = sp.zeros(len(x))
        for i in xrange(n+1):
            z += c[i] * x**(n-i)
        z += c[-1] * xii
        dmagseg1 = mag[seg1] - z + np.median(mag)
        
        seg2 = data['cadence'] > 624
        x = data['cadence'][seg2]
        xi = sp.zeros(len(x)) + 1
        xii = sp.zeros(len(x)) + focus[seg2]
        amatrix = (x**n)[:, np.newaxis]
        for i in xrange(1, n):
            amatrix = np.hstack((amatrix, (x**(n-i))[:, np.newaxis]))
        amatrix = np.c_[amatrix, xi[:, np.newaxis]]
        amatrix = np.c_[amatrix, xii[:, np.newaxis]]
        c, resid, rank, sigma = sp.linalg.lstsq(amatrix, noflux[seg2])
        z = sp.zeros(len(x))
        for i in xrange(n+1):
            z += c[i] * x**(n-i)
        z += c[-1] * xii
        dmagseg2 = mag[seg2] - z + np.median(mag)
        dmag[seg1] = dmagseg1
        dmag[seg2] = dmagseg2
        return dmag


def lspolyordern(x, y, n):
    x = sp.array(x)
    y = sp.array(y)
    xi = sp.zeros(len(x)) + 1
    amatrix = (x**n)[:, np.newaxis]
    for i in range(1, n):
        amatrix = np.hstack((amatrix, (x**(n-i))[:, np.newaxis]))

    amatrix = np.c_[amatrix, xi[:, np.newaxis]]
    c, resid, rank, sigma = sp.linalg.lstsq(amatrix, y)
    z = sp.zeros(len(x))
    for i in range(n+1):
        z += c[i]*x**(n-i)
    return [z, c]


def matrixgen(time, n, timespan):
    # generate least square fitting matrix with n cos filter,
    # t is the total time span, formulism refer to Huang and Bacos (2012) eq [1]'''
    tn = len(time)
    a = np.rollaxis(np.sin(np.array(np.r_['c', 0:n] * time[np.arange(tn)])
                           * math.pi / timespan), 1, 0)
    return a


def flatfunc(c, ctime, timespan):
    n = len(c) - 1
    b = np.rollaxis(np.sin(np.array(np.r_['c', 0:n] * ctime) * math.pi / timespan),
                    1, 0)
    rflux = np.dot(b, c[0:n])
    rflux -= np.mean(rflux)
    return rflux


def tranfunc(c, ctime, intran, timespan):
    n = len(c)-1
    b = np.rollaxis(np.sin(np.array(np.r_['c', 0: n] * ctime) * math.pi / timespan),
                    1, 0)
    rflux = np.dot(b, c[0: n])
    try:
        rflux[intran] += c[n]
    except TypeError:
        print c
        raise
    rflux -= np.mean(rflux)
    return rflux


def lsfitrecon(otime, oflux, intran, n=30, noplot=True, dipguess=0.0001):
    length = len(oflux)
    ctime = np.array(otime)
    epoch0 = min(ctime)
    timespan = max(ctime)-min(ctime)
    cflux = np.array(oflux)-np.mean(oflux)

    def e(c, time, index, flux, timespan):
        return tranfunc(c, time, index, timespan) - flux

    c0 = list(np.zeros(n))
    c0.append(dipguess)
    c, success = leastsq(e, c0, args=(ctime - epoch0, intran, cflux, timespan), maxfev=10000)
    b = np.rollaxis(np.sin(np.array(np.r_['c', 0:n]
                                    * (ctime[np.arange(length)]-E0))
                           * math.pi/timespan), 1, 0)
    rflux = np.dot(b, c[0:n])
    tran = sp.zeros(len(ctime))
    tran[intran] += c[n]
    tran -= np.mean(tran)
    # print c[n],len(tran[intran])
    if not noplot:

        import matplotlib
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(ctime, cflux, 'x', ctime, rflux, '.', ctime, tran, '+')
        plt.show()
        plt.close(fig)
    # rl = len(rflux)
    # nn = len(otime)
    dflux = oflux - rflux + np.mean(oflux)
    return dflux
