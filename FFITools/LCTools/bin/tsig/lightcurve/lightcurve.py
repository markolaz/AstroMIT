#
# Copyright (C) 2015 - Zach Berta-Thompson <zkbt@mit.edu> (MIT License)
#               2016 - Chelsea Huang (GPLv3 License)
#               2016 - Massachusetts Institute of Technology (MIT)
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
Base class for light curves and light curve generation classes.
"""

import ConfigParser
import logging
import random
import numpy as np
import math
import os
from tsig.util import constants
from math import pi
#from dataio import readcolumn

logger = logging.getLogger(__name__)

def create_random_lightcurve(options=('trapezoid', 'sin'),
                             fraction_with_extreme=0.01,
                             fraction_with_rotation=None,
                             fraction_with_trapezoid=None,
                             fraction_with_custom=0.0):
    """
    Randomly choose a type of light curve.

    options=(trapezoid | sin)
        a list of the kinds of a variability to choose from)

    fraction_with_extreme=0.01
        what fraction should be drawn from a separate "extreme" population
        of light curves?

    fraction_with_trapezoid=None
        what fraction should get trapezoid? if None, defaults to 18%,
        from Kepler

    fraction_with_rotation=None
        what fraction should get rotation? if None, defaults to 26%,
        from Kepler

    when injecting, a star gets (in this order):
        - a chance to be extreme
        - a chance to be a trapezoid
        - a chance to be a rotator
        """

    # first of all, give preference to try to be extreme
    if np.random.uniform(0, 1) < fraction_with_extreme:
        return create_cartoon_lightcurve(options=options, extreme=True)
    else:
        # set the fraction that get rotation to that from Kepler
        if fraction_with_rotation is None:
            fraction_with_rotation = 34030.0 / 133030.0
        # set the fraction that get transits to that from Kepler
        if fraction_with_trapezoid is None:
            fraction_with_trapezoid = 20152 / 112001.0
        # give trapezoids preference over sin curves...
        if 'trapezoid' in options:
            if np.random.uniform(0, 1) < fraction_with_trapezoid:
                return KeplerTCELightCurve()
        # ... otherwise try to include a sine curve
        if 'sin' in options:
            if np.random.uniform(0, 1) < fraction_with_rotation:
                return McQuillanLightCurve()
    # if nothing else, make the light curve a constant
    return ConstantLightCurve()


def create_cartoon_lightcurve(options=('trapezoid', 'sin'),
                              extreme=False,
                              prng=np.random):
    """Generate a random lightcurve, from a cartoonish population"""

    name = prng.choice(options)
    if extreme or name == 'sin':
        if extreme:
            p = [0.1, 30.0]
            a = [0.1, 1]
        else:
            p = [0.1, 30.0]
            a = [0.0001, 0.02]

        P = 10 ** prng.uniform(*np.log10(p))
        E = prng.uniform(0, P)
        A = 10 ** prng.uniform(*np.log10(a))

        return SinusoidLightCurve(P=P, E=E, A=A)

    elif name == 'trapezoid':
        if extreme:
            p = [0.1, 30.0]
            d = [0.1, 1]
        else:
            p = [0.1, 30.0]
            d = [0.0001, 0.01]
        P = 10 ** prng.uniform(*np.log10(p))
        E = prng.uniform(0, P)

        mass = prng.uniform(0.1, 1.5)
        radius = mass
        stellar_density = 3 * mass * constants.Msun / (4 * np.pi * (radius * constants.Rsun) ** 3)
        rsovera = (3 * np.pi / constants.G / (P * constants.day) ** 2 / stellar_density) ** (1.0 / 3.0)
        T14 = rsovera * P / np.pi
        T23 = prng.uniform(0, T14)
        D = 10 ** prng.uniform(*np.log10(d))

        return TrapezoidLightCurve(P=P, E=E, D=D, T23=T23, T14=T14)

    else:
        raise ValueError("No option specified.  Must be 'trapezoid' or 'sin'.")

class classproperty(object):
    def __init__(self, f):
        self.f = f

    def __get__(self, obj, owner):
        return self.f(owner)


class LightCurve(dict):
    """
    The LightCurve class defines the basics of a light curve object, which can
    injected into TESS simulations. It handles basic functionality, like
    (importantly), integrating a light curve over a finite exposure time.
    """

    def __repr__(self):
        """string representation of the light curve object"""
        return '<{0}>'.format(self.code)

    @property
    def code(self):
        """returns a string describing the lightcurve"""
        t = ''
        for k, v in self.iteritems():
            try:
                t += '{k}={v:0.3f},'.format(k=k, v=v)
            except TypeError:
                t += '{k}={v},'.format(k=k, v=v)
        return '{0}({1})'.format(self.__class__.__name__, t[:-1])

    def get_flux(self, t):
        raise NotImplementedError()

    def get_mag(self, t):
        """
        Return an array of magnitudes whose values are relative to a baseline.
        """
        raise NotImplementedError

    def make_plot(self, tmin=0, tmax=27.4,
                  cadence=30.0 / 60.0 / 24.0, offset=0, raw=False, ax=None):
        """make a plot of a light curve, to show what it looks like"""
        import matplotlib.pyplot as plt
        if ax is None:
            plt.figure(str(self), figsize=(8, 3))
        else:
            ax.set_title(str(self))
            plt.sca(ax)
        t = np.arange(tmin, tmax, cadence)
        y = self.get_mag(t)
        if raw:
            plt.plot(t, y + offset, alpha=0.25, linewidth=1, color='royalblue')
        plt.plot(t, self.integrated(t) + offset, alpha=0.5,
                 linewidth=1, color='darkorange')
        plt.xlim(tmin, tmax)
        # plt.ylim(np.max(y)+0.01, np.min(y)-0.01)
        plt.xlabel('Time (days)')
        plt.ylabel('Flux (mag)')

        # Return the current working figure and NOT the plt class, this is so
        # we don't end up with concurancy issues if other figures are created
        # and our calling code doesn't understand why their graph got deleted.
        return plt.gcf()

    def chelseaplot(self, show=True, label='rlc'):
        # only works for the sub
        import matplotlib
        from matplotlib import pyplot as plt
        self.check_data('jd')
        self.check_data(label)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if self.data['jd'][0] > 2454000:
            ax.plot(self.data['jd'] - 2454000, self.data[label], 'k.')
            ax.set_xlabel("BJD-2454000")
        else:
            ax.plot(self.data['jd'], self.data[label], 'k.')
            ax.set_xlabel("BJD")
        
        ax.set_ylabel(self.labels[label])
        yformatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(yformatter)
        ax.grid(True)
        if show:
            plt.show()
        else:
            pngfile = os.path.splitext(self.name)[0] + '.png'
            plt.savefig(pngfile)

    def integrated(self, t, exptime=30.0 / 60.0 / 24.0, resolution=100):
        """
        Integrate the flux over a finite exposure time.

          t - one or more time stamps
          exptime - duration in days
          resolution - 

        Returns an array of values with shape equal to the input t.

        For proper averaging, the values returned by the lightcurve model
        must be defined in magnitudes.
        """

        # FIXME: why is the result of integration a vector?
        
        # convert a single time point to an array with one element
        try:
            t.shape
        except AttributeError:
            if isinstance(t, list) or isinstance(t, tuple):
                t = np.array(t)
            else:
                t = np.array([t])

        # constant light curve is a degenerate case
        if self.__class__.__name__ == 'ConstantLightCurve':
            return self.get_mag(np.array(t))

        # create a high-resolution subsampled timeseries
        nudges = np.linspace(-exptime / 2.0, exptime / 2.0, resolution)
        subsampled = t.reshape(1, t.shape[0]) + nudges.reshape(nudges.shape[0], 1)

        # make sure the average is photon-weighted (not magnitude-weighted)
        flux = 10 ** (-0.4 * self.get_mag(subsampled))
        mag = -2.5 * np.log10(flux.mean(0))
        assert (mag.shape == t.shape)
        return mag

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

    def __getstate__(self):
        d = self.__dict__.copy()
        if 'logger' in d.keys():
            d['logger'] = d['logger'].name
        return d

    def __setstate__(self, d):
        if 'logger' in d.keys():
            d['logger'] = logging.getLogger(d['logger'])
        self.__dict__.update(d)


class ConstantLightCurve(LightCurve):
    """Light curve whose brightness is constant."""

    def __init__(self, offset=0):
        super(ConstantLightCurve, self).__init__(offset=offset)

    def get_mag(self, t):
        return np.zeros_like(t) + self['offset']


class SinusoidLightCurve(LightCurve):
    """
    Light curve whose brightness is defined by a single sinusoid defined by
    period, phase, and amplitude.
    """

    def __init__(self, P=pi, E=0.0, A=0.1):
        super(SinusoidLightCurve, self).__init__(P=P, E=E, A=A)

    def get_mag(self, t):
        return self['A'] * np.sin(2 * np.pi * (t - self['E']) / self['P'])


class TrapezoidLightCurve(LightCurve):
    """
    Light curve whose brightness is defined by a simplified eclipse model:
        P = period
        E = epoch of one eclipse center
        D = depth
        T14 = duration between contacts 1-4
        T23 = duration between contacts 2-3
    """

    def __init__(self, P=pi, E=0.0, D=0.01, T23=0.1, T14=0.1):
        super(TrapezoidLightCurve, self).__init__(
            P=P, E=E, D=D, T23=T23, T14=T14)

    def closest_transit(self, t):
        return np.round((t - self['E']) / self['P']) * self['P'] + self['E']

    def timefrommidtransit(self, t):
        return t - self.closest_transit(t)

    def get_mag(self, t):
        """model is returned in magnitudes, relative to a baseline level"""
        flux = np.zeros_like(t)
        dt = np.abs(self.timefrommidtransit(t))
        start, finish = self['T23'] / 2.0, self['T14'] / 2.0
        depth = self['D']

        i = dt <= start
        flux[i] = depth
        i = (dt <= finish) * (dt > start)
        flux[i] = (finish - dt[i]) / (finish - start) * depth

        return flux


class McQuillanLightCurve(SinusoidLightCurve):
    """Random sine curve, drawn from McQuillan rotation periods"""

    def __init__(self, prng=np.random):
        row = prng.choice(McQuillanLightCurve.get_rotation_table())
        super(McQuillanLightCurve, self).__init__(
            P=row['PRot'],
            E=prng.uniform(0, row['PRot']),
            A=row['Rper'] / 1.0e6)

    @staticmethod
    def get_rotation_table():
        """Lazy load the McQuillan rotation periods table.  Once loaded, the
        table is retained as a static data member of the LightCurve class."""
        if not hasattr(McQuillanLightCurve, 'rotation_table'):
            import astropy.io.ascii as ascii
            import pkgutil
            logger.debug("Reading McQuillan rotation table")
            McQuillanLightCurve.rotation_table = ascii.read(
                pkgutil.get_data(__name__, 'rotation_McQuillan.txt'))
        return McQuillanLightCurve.rotation_table


class KeplerTCELightCurve(TrapezoidLightCurve):
    """Random transit from the Kepler TCE list"""

    def __init__(self, prng=np.random):
        row = prng.choice(KeplerTCELightCurve.get_transit_table())
        T14 = row['tce_duration'] / 24.0
        super(KeplerTCELightCurve, self).__init__(
            P=row['tce_period'],
            E=row['tce_time0bk'] + 2545833.0,
            T14=T14,
            T23=max(T14 - 2 * row['tce_ingress'] / 24.0, 0),
            D=row['tce_depth'] / 1.0e6)

    @staticmethod
    def get_transit_table():
        """Lazy load the Kepler TCE table.  Once loaded, the table is retained
        as a static data member of the LightCurve class."""
        # from Seader et al.
        #   "The search includes a total of $198,646$ targets, of which
        #   $112,001$ were observed in every quarter and $86,645$ were
        #   observed in a subset of the 17 quarters."
        if not hasattr(KeplerTCELightCurve, 'transit_table'):
            import astropy.io.ascii as ascii
            import pkgutil
            logger.debug("Reading Kepler TCE table")
            KeplerTCELightCurve.transit_table = ascii.read(
                pkgutil.get_data(__name__, 'keplerTCE_DR24.txt'))
        return KeplerTCELightCurve.transit_table


class RandomLightCurve(LightCurve):
    """Light curve with randomly generated flux values using a normal
    distribution"""

    def __init__(self, sigma=0.001, depth=0.01, CD=0.02083333333,
                 E=2, H=0.4, P=10, m0=1):
        super(RandomLightCurve, self).__init__(sigma=sigma, depth=depth,
                                               CD=CD, E=E, H=H, P=P, m0=m0)

    def get_flux(self, t):
        # Random light curve generated using technique in andras' se-mklc.sh
        # FIXME: make this handle vector t
        v = random.gauss(0, 1)
        phase = (t - self['E']) / self['P']
        phase -= int(phase)
        while phase < -0.5:
            phase += 1
        while 0.5 < phase:
            phase -= 1
        dt = phase * self['P']
        mag = self['m0']
        if -self['H'] <= dt <= self['H']:
            mag += self['depth']
        return mag + v * self['sigma']

    def get_mag(self, t):
        # FIXME: this does not match logic of get_flux
        flux = np.zeros_like(t)
        flux += self['m0']
        flux += np.random.random(flux.shape) * self['sigma']
        return flux


class FilebasedLightCurve(LightCurve):
    """Light curve based on user input file"""

    Z0 = 20.17 #tess

    def __init__(self, name=''):
        self.name = name
        self.data = {'jd': [], # Julian Date
                     'rlc': [], # Raw light curve data
                     'ltflc': [], #Processed light curve data
                     'x': [], # Pixel x-position (used for jitter)
                     'y': [], # Pixel y-position (used for jitter)
                     'bg': [], # Star background measurements
                     'bgerr': [], # Star background error
                     'cadence': []} # Frame number
        self.labels = {'jd': 'BJD', 'rlc': 'Raw LC (Magnitude)',
                       'ltflc': 'DLC (Magnitude)', 'x': 'X (Pixel)',
                       'y': 'Y (Pixel)', 'bg': 'Background',
                       'bgerr': 'Background Err', 'cadence': 'Cadence'}
        self.logger = logging.getLogger(LightCurve.__name__)
        self.logger.addHandler(logging.NullHandler())
        return
    
    @staticmethod
    def create_lightcurve(infile,cfgfile='example.cfg', useridentifier=None):
        fin = open(infile,mode='r')
        identifier = fin.read(8)
        fin.close() 
        for c in [FITSLightCurve, HDFLightCurve, ASCIILightCurve]:
            if c.istype(identifier):
                return c(name=infile,cfgfile=cfgfile)
        if useridentifier:
            for c in [FITSLightCurve, HDFLightCurve, ASCIILightCurve]:
                if c.istype(useridentifier):
                    return c(name=infile,cfgfile=cfgfile)

        raise Exception, 'File format not recognized for input file %s' % infile 

    def load_from_file(self, label='all'):
        """
        Load data from a file.  Allow user to specify load all or specific
        columns.
        """
        raise NotImplementedError

    def write_to_file(self, outfile):
        # Output lightcurve to a file as designed.
        raise NotImplementedError

    def _gen_jd(self):
        self.check_data(label='cadence')
        self.data['jd'] = self.data['cadence'] / 48. + 2457827.0

    def check_data(self, label='rlc'):
        if len(self.data[label]) == 0:
            try:
                self.load_from_file(label=label)
            except (IndexError, KeyError):
                raise IndexError, "Column %s is required, but cannot found in file %s" % (label, self.name)

    def get_flux(self, t):
        try:
            return 10**((FilebasedLightCurve.Z0 - self.get_mag(t))/2.5)/10**((FilebasedLightCurve.Z0 - sum(self.data['rlc'])/len(self.data['rlc']))/2.5)
        except TypeError:
            return None

    def get_mag(self, t):
        time, mag = self.data['jd'], self.data['rlc']
        minutetime = [(time[x]-time[0])*24*60 for x in xrange(len(time))]
        try:
            return mag[minutetime.index(t)]
        except ValueError:
            try:
                lowerlimit = min(minutetime, key=lambda x:abs(x-t))
                upperlimit = minutetime[minutetime.index(lowerlimit)+1]
                t1, t2 = minutetime[minutetime.index(lowerlimit)], \
                         minutetime[minutetime.index(upperlimit)]
                m1, m2 = mag[minutetime.index(lowerlimit)], \
                         mag[minutetime.index(upperlimit)]
                return m1 + (t - t1)*(m2-m1)/(t2-t1)
            except IndexError:
                return None

class ASCIILightCurve(FilebasedLightCurve):
    def __init__(self, cfgfile='example.cfg', name=''):
        super(ASCIILightCurve, self).__init__(name)
        self.logger = logging.getLogger(ASCIILightCurve.__name__)
        self.logger.addHandler(logging.NullHandler())
        p = ConfigParser.RawConfigParser()
        p.read(cfgfile)
    
        try: 
            self.cols = eval(p.get('LC', 'lcformat'))
            self.logger.debug("lcformat %s", str(self.cols))
        except ConfigParser.NoOptionError:
            self.cols = {'jd': 2, 'rlc': 8, 'ltflc': 12,
                        'x': 3, 'y': 4, "bg": 12, "bgerr": 12, 'cadence': 2}
            self.logger.warning("Use default LC format, type is %s", str(self.cols))  
        except ConfigParser.NoSectionError:
            raise ConfigParser.NoSectionError, "No LC section in configuration file %s" % cfgfile
        
        self.outext = ''
    
    @staticmethod 
    def istype(identifier):

        if all(ord(c) <128 for c in identifier): 
            return True
        return False

    def load_from_file(self, label='all'):
        with open(self.name, mode='r') as fin:
            headerline = fin.readline()
        
        if headerline.startswith('#'):
            newdict = dict(zip(headerline.split()[1:],np.arange(len(headerline.split())-1)+1))
            self.cols.update(newdict)
            self.logger.debug("Using columns from file: %s", str(self.cols))
        else: 
            self.logger.debug("File has no header line, using default format or format from configuration file: %s", str(self.cols))

        if label == 'all':
            self.readcolumn(self.data['jd'], int(self.cols['jd']), self.name)
            self.data['jd'] = np.array(self.data['jd'])
            self.readcolumn(self.data['rlc'], int(self.cols['rlc']), self.name)
            self.data['rlc'] = np.array(self.data['rlc'])
            try:
                self.readcolumn(self.data['ltflc'], int(self.cols['ltflc']), self.name)
                self.data['ltflc'] = np.array(self.data['ltflc'])
            except (IndexError, KeyError): 
                self.logger.warning("ltflc for %s not created yet, you might want to detrend first", self.name)
            try:
                self.readcolumn(self.data['x'], int(self.cols['x']), self.name)
                self.data['x'] = np.array(self.data['x'])
                self.readcolumn(self.data['y'], int(self.cols['y']), self.name)
                self.data['y'] = np.array(self.data['y'])
            except (IndexError, KeyError): 
                self.logger.warning("coordinate columns x, "
                                    "y for %s not created yet, you might want "
                                    "to detrend first", self.name)
            try:
                self.readcolumn(self.data['bg'], int(self.cols['bg']), self.name)
                self.data['bg'] = np.array(self.data['bg'])
                self.readcolumn(self.data['bgerr'], int(self.cols['bgerr']), self.name)
                self.data['bgerr'] = np.array(self.data['bgerr'])
            except (IndexError, KeyError): 
                self.logger.warning("Warning: background columns bg, bgerr "
                                    "for %s not created yet, you might want "
                                    "to detrend first", self.name)
        else:
        
            self.readcolumn(self.data[label], int(self.cols[label]), self.name)
            self.data[label] = np.array(self.data[label])


    def write_to_file(self, outfile, replace=True, outheader="# jd rlc"):
        # Output lightcurve to a file as designed.
        if os.path.exists(outfile) and (not replace):
            self.logger.warning("write aborted: file exists at %s", outfile)
            return
        self.logger.info("outfile header is %s", outheader)
        with open(outfile, mode='w') as fout:
            fout.write(outheader)
            fout.write('\n')        
            for i in xrange(len(self.data['jd'])):
                if np.isnan(self.data['rlc'][i]):
                    continue
                for label in outheader.split()[1:]:
                    fout.write("%15.7f " % (self.data[label][i]))
                fout.write("\n")

    def readcolumn(self,var,col,infile,datformat='float',div=None,fill=False):
        with open(infile,mode='r') as fin:
            data=fin.readline().split(div)
            while data:           
                if not data[0].startswith('#'):
                    try:
                        if datformat=='float':
                            var.append(float(data[col-1]))
                        else:
                            var.append(data[col-1])
                    except (ValueError, IndexError):
                        if fill:
                            if datformat=='float':
                                var.append(np.nan)
                            else:
                                var.append(None)
                        else:
                            raise
                data=fin.readline().split(div)


class HDFLightCurve(FilebasedLightCurve):
    # TBD. need to merge from old branch
    IDENTIFIER = "\211HDF\r\n\032\n"
    def __init__(self,name=''):
        super(HDFLightCurve, self).__init__(name)
        self.logger = logging.getLogger(HDFLightCurve.__name__)
        self.logger.addHandler(logging.NullHandler())
        return
    
    @staticmethod
    def istype(identifier):
        if identifier == HDFLightCurve.IDENTIFIER:
            return True
        return  False
    
    def load_from_file(self, label='all'):
        raise NotImplementedError

    def write_to_file(self, outfile, replace=True, outheader="# jd rlc"):
        raise NotImplementedError


class FITSLightCurve(FilebasedLightCurve):
    IDENTIFIER = "SIMPLE  "
    BJD = 2454000
    #Z0 = 24.2415 #Kepler
    Z0 = 20.17 #tess
    def __init__(self, name='', cfgfile='example.cfg'):
        super(FITSLightCurve, self).__init__(name)
        self.logger = logging.getLogger(FITSLightCurve.__name__)
        self.logger.addHandler(logging.NullHandler())

    @staticmethod
    def istype(identifier):
        if identifier == FITSLightCurve.IDENTIFIER:
            return True
        return False
    def load_from_file(self, label='all'):
        # load data from a file, allow user to specify load all or specific columns.
        import astropy
        import astropy.io.fits as fits 
        with fits.open(self.name) as hdulist:
            tbdata = hdulist[1].data
            
            BJDREFI = hdulist[1].header['BJDREFI']
            BJDREFF = hdulist[1].header['BJDREFF']
     
            self.data['jd'] = np.array(tbdata.field('Time')) + BJDREFI + BJDREFF - FITSLightCurve.BJD
            self.logger.debug("lenth of data is %d", len(self.data['jd']))
            self.data['rlc'] = FITSLightCurve.Z0 - 2.5*np.log10(np.array(tbdata.field('PDCSAP_FLUX')))
            try:
                TESSmag=hdulist[0].header['TESSMAG']
                self.data['rlc']-=(np.nanmedian(self.data['rlc'])-TESSmag)
            except KeyError:
                pass
            self.data['x'] = np.array(tbdata.field('MOM_CENTR1'))
            self.data['y'] = np.array(tbdata.field('MOM_CENTR2')) 
            self.data['cadence'] = np.array(tbdata.field('CADENCENO')) 

    def write_to_file(self, outfile, replace=True, outheader="# jd rlc"):
        # Output lightcurve to a file as designed.
        if os.path.exists(outfile) and (not replace):
            self.logger.warning("write aborted: file exists at %s", outfile)
            return
        self.logger.info("outfile header is %s", outheader)
        with open(outfile, mode='w') as fout:
            fout.write(outheader)
            fout.write('\n')        
            for i in xrange(len(self.data['jd'])):
                if np.isnan(self.data['rlc'][i]):
                    continue
                for label in outheader.split()[1:]:
                    fout.write("%15.7f " % (self.data[label][i]))
                fout.write("\n")
