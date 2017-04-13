#!/usr/bin/env python
import ConfigParser
import optparse
import logging
import sys
import os
import numpy as np

import logging
import sys

try: 
    import Tkinter
    mpl_backend = 'TkAgg'
except ImportError:
    mpl_backend = 'Agg'

import matplotlib
matplotlib.use(mpl_backend)
from matplotlib import pyplot as plt
from matplotlib import rc


from lctools.util.util import getext, mad 
from lctools.util.configurable import ConfigurableObject
from lctools.util.dataio import readtableline
from lctools.transit import Transit
# from lctools.transit import calc_epoch, get_fold_fit, mad, transit_fit, modeltransit, occultquad 
# create logger
logger = logging.getLogger(__name__)

class Summary(object):
    def __init__(self, lc, blsanal, **kwargs):
        super(Summary, self).__init__()
        
        self.figure = plt.figure(figsize=[16,9])
        self.setup_plot()
        self.required_cols = ['jd', 'ltflc', 'x', 'y']
        self.check_data(lc)
                
        index = np.isnan(lc.data['ltflc'])
        self.jd = lc.data['jd'][~index]
        self.flux = 10.**(-0.4*(lc.data['ltflc'][~index]-np.nanmedian(lc.data['ltflc'])))
        
        self.x = lc.data['x'][~index]
        self.y = lc.data['y'][~index]
        
        candparam = readtableline(blsanal)
        self.transit = Transit(name=os.path.splitext(os.path.basename(lc.name))[0].split('_')[0], P=candparam['Period'], q=candparam['Qtran'], qg = candparam['Qingress'], depth = candparam['Depth'], epoch = candparam['Tc'], SN=candparam['SN'], SNR=candparam['SignaltoPinknoise'])
        
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
 
    def check_data(self, lightcurve):
        for key in self.required_cols:
            lightcurve.check_data(key)
        return 

    def setup_global_parameters(self):
        rc('text', usetex=True)
        rc('font', family='serif')
        rc('font', serif='cm')
        rc('font', size=18, weight='black', stretch='extra-expanded', style='oblique')
        rc('axes', linewidth=2)

    def setup_plot(self):
        self.setup_global_parameters() 
    
        # detrended full light curve plot (top)
        self.ax1 = self.figure.add_subplot(331)
        pos1 = self.ax1.get_position()
        pos2 = [pos1.x0, pos1.y0+0.03, pos1.width*2.0+0.03, pos1.height] 
        self.ax1.set_position(pos2)

        self.ax2 = self.figure.add_subplot(334)

        #Secondary Eclipse
        self.ax3 = self.figure.add_subplot(335)
        pos1 = self.ax3.get_position()
        pos2 = [pos1.x0+0.03, pos1.y0, pos1.width, pos1.height]
        self.ax3.set_position(pos2)

        #Even Transit
        self.ax4 = self.figure.add_subplot(337)
        pos1 = self.ax4.get_position()
        pos2 = [pos1.x0-0.02, pos1.y0-0.03, pos1.width, pos1.height] 
        self.ax4.set_position(pos2)
        
        #Odd Transit
        self.ax5 = self.figure.add_subplot(338)
        pos1 = self.ax5.get_position()
        pos2 = [pos1.x0-0.02, pos1.y0-0.03, pos1.width, pos1.height] 
        self.ax5.set_position(pos2)

        #Centroid X Plot
        self.ax6 = self.figure.add_subplot(339)
        pos1 = self.ax6.get_position()
        pos2 = [pos1.x0+0.03, pos1.y0-0.03, pos1.width, pos1.height] 
        self.ax6.set_position(pos2)
        
        #Centroid Y Plot
        self.ax7 = self.ax6.twinx()
        self.ax7.set_position(pos2)
    
    def gen_summary(self, plottofile=False):
        
        self.plot_figure1()
        self.plot_figure2()
        self.plot_figure3()
        self.plot_figure4to5()
        self.plot_figure6to7()

        
        if plottofile:
            figfile = os.path.splitext(self.transit.name)[0]+'.png'
            plt.savefig(figfile, dpi=150)
        else:
            plt.show()

    def plot_figure1(self):
        self.ax1.set_xlabel("BJD")
        self.ax1.set_ylabel("Relative Flux (ppm)")
        self.ax1.annotate("Detrended Full LC", [0.3, 0.95], xycoords='figure fraction',  fontweight='black', color='#be7013', fontsize=18, fontstyle='oblique')
   
        self.ax1.plot(self.jd, (self.flux-1)/1.e-6, lw=0, marker='.', color='k')
        sigma = mad(self.flux)
        self.ax1.set_ylim((min(self.flux)-1.)/1.e-6, (7*sigma)/1.e-6)  # inverted y-axis
        self.ax1.set_xlim(min(self.jd), max(self.jd))  # inverted y-axis
        self.ax1.get_yaxis().get_major_formatter().set_useOffset(False)

        # end = int(np.floor((self.jd[-1] - self.Tc) / self.P) + 1)
        # start = int(np.floor((self.jd[0] - self.Tc) / self.P))
        # logger.debug("end=%d, start=%d, t[0] = %f, t[-1] = %f, p=%f", end, start, self.jd[0], self.jd[-1], self.P)


        self.transit.calc_epoch(self.jd)

        for midpt in self.transit.midpts:
            self.ax1.vlines(midpt, (min(self.flux)-1.)/1.e-6, (max(self.flux)-1.)/1.e-6 , color='r', linestyle='--')



        #annotate basic information of the candidate
        self.ax1.annotate("Planet: $%s$" % self.transit.name, [0.7, 0.93], xycoords='figure fraction',  fontweight='black', color='k', fontsize=25, fontstyle='oblique')
        self.ax1.annotate("$P = %.3f {\ \mathrm {Day}}$" % self.transit.P, [0.7, 0.88], xycoords='figure fraction',  fontweight='black', color='k', fontsize=20, fontstyle='oblique')
        self.ax1.annotate("$T_0 = %.3f {\ \mathrm {BJD}}$" % self.transit.epoch, [0.7, 0.83], xycoords='figure fraction',  fontweight='black', color='k', fontsize=20, fontstyle='oblique')
        self.ax1.annotate("$Rp/Rstar \sim %.3f$" % np.sqrt(self.transit.depth), [0.7, 0.78], xycoords='figure fraction',  fontweight='black', color='k', fontsize=20, fontstyle='oblique')
        self.ax1.annotate("$T_{dur}/P \sim %.3f$" % self.transit.q, [0.7, 0.73], xycoords='figure fraction',  fontweight='black', color='k', fontsize=20, fontstyle='oblique')
        self.ax1.annotate("$T_{12}/T_{14} \sim %.3f$" % self.transit.qgress, [0.7, 0.68], xycoords='figure fraction',  fontweight='black', color='k', fontsize=20, fontstyle='oblique')
        self.ax1.annotate("$SN_{\mathrm {BLS}} \sim %.3f$" % self.transit.SN, [0.7, 0.63], xycoords='figure fraction',  fontweight='black', color='k', fontsize=20, fontstyle='oblique')
        self.ax1.annotate("$SNR \sim %.3f$" % self.transit.SNR, [0.7, 0.58], xycoords='figure fraction',  fontweight='black', color='k', fontsize=20, fontstyle='oblique')
        self.ax1.annotate("Star: KICxxxxx ", [0.7, 0.53], xycoords='figure fraction',  fontweight='black', color='k', fontsize=25, fontstyle='oblique')
        self.ax1.annotate("$R_* \sim R_{\odot}$", [0.7, 0.48], xycoords='figure fraction',  fontweight='black', color='k', fontsize=20, fontstyle='oblique')
        self.ax1.annotate("$M_* \sim M_{\odot}$", [0.8, 0.48], xycoords='figure fraction',  fontweight='black', color='k', fontsize=20, fontstyle='oblique')
        self.ax1.annotate("$T_{\mathrm {eff}} \sim K$", [0.7, 0.43], xycoords='figure fraction',  fontweight='black', color='k', fontsize=20, fontstyle='oblique')
        self.ax1.annotate("$logg \sim $", [0.8, 0.43], xycoords='figure fraction',  fontweight='black', color='k', fontsize=20, fontstyle='oblique')
        self.ax1.annotate("${\mathrm {RA}} =  $", [0.7, 0.38], xycoords='figure fraction',  fontweight='black', color='k', fontsize=20, fontstyle='oblique')
        self.ax1.annotate("${\mathrm {DEC}} =  $", [0.8, 0.38], xycoords='figure fraction',  fontweight='black', color='k', fontsize=20, fontstyle='oblique')
        self.ax1.annotate("${\mathrm {KepMag}} =  $", [0.7, 0.33], xycoords='figure fraction',  fontweight='black', color='k', fontsize=20, fontstyle='oblique')
        self.ax1.annotate("$J-K =  $", [0.8, 0.33], xycoords='figure fraction',  fontweight='black', color='k', fontsize=20, fontstyle='oblique')
        self.ax1.annotate("${\mathrm {DEC}} =  $", [0.8, 0.38], xycoords='figure fraction',  fontweight='black', color='k', fontsize=20, fontstyle='oblique')
        
    def plot_figure2(self):
        self.ax2.set_xlabel("Time from T0 (hours)")
        self.ax2.set_ylabel("Relative Flux (ppm)")
        self.ax2.annotate("Primary Transit", [0.14, 0.4], xycoords='figure fraction',  fontweight='black', color='#be7013', fontsize=18, fontstyle='oblique')
        dt_tra, f_tra, binned_t, binned_flux = self.transit.get_primary(self.jd, self.flux) 
        sigma = mad(f_tra)
        self.ax2.plot(dt_tra*24, f_tra/1.e-6, lw=0, marker='.', color='0.7')
        self.ax2.plot(binned_t*24, binned_flux/1e-6, marker='.', lw=0, color='#1c79e8', markersize=20)
        self.ax2.axvline(x=0, ls='--', color='k')
        self.ax2.set_ylim((min(f_tra)-np.median(f_tra))/1.e-6, sigma*6/1.e-6)
        self.ax2.get_yaxis().get_major_formatter().set_useOffset(False)
        
    def plot_figure3(self):
        self.ax3.set_xlabel("Time from T0 (hours)")
        self.ax3.set_ylabel("Relative Flux (ppm)")
        # ax3.set_title('Occultation', color='k', fontsize=18)
        self.ax3.annotate("Occultation", [0.44, 0.4], xycoords='figure fraction',  fontweight='black', color= '#be7013', fontsize=18, fontstyle='oblique')
        t_occ, f_occ, t_bin, f_bin = self.transit.get_occultation(self.jd, self.flux) 
        sigma = np.std(f_occ)

        self.ax3.plot(t_occ*24., f_occ/1.e-6, lw=0, marker='.', color='0.75')
        self.ax3.plot(t_bin*24., f_bin/1.e-6, lw=2, marker='.', color='#1c79e8', markersize=20)
        self.ax3.set_ylim(-3*sigma/1.e-6, 2*sigma/1.e-6)
        self.ax3.get_yaxis().get_major_formatter().set_useOffset(False)

    def plot_figure4to5(self):
        self.ax4.set_xlabel("Time from T0 (hours)")
        self.ax4.set_ylabel("Relative Flux (ppm)")
        self.ax4.annotate("Even Transit", [0.1, 0.1], xycoords='figure fraction',  fontweight='black', color= '#be7013', fontsize=18, fontstyle='oblique')
    
       
        self.ax5.set_xlabel("Time from T0 (hours)")
        self.ax5.get_yaxis().set_ticks([])
        self.ax5.annotate("Odd Transit", [0.38, 0.1], xycoords='figure fraction',  fontweight='black', color= '#be7013', fontsize=18, fontstyle='oblique')
        
        dt_odd, f_odd, dt_odd_bin, f_odd_bin, tarr_odd, oddmod, dt_even, f_even, dt_even_bin, f_even_bin, tarr_even, evenmod, sigma = self.transit.get_odd_even(self.jd, self.flux)
        odd_depth = np.min(oddmod) 
        even_depth = np.min(evenmod) 
        difference = np.abs(even_depth-odd_depth)/sigma
        self.ax4.plot(dt_odd * 24., f_odd/1.e-6, lw=0, marker='.', color='0.7')
        self.ax4.plot(24*dt_odd_bin, f_odd_bin/1e-6, marker='.', color='#1c79e8', markersize=20, lw=0)
        self.ax4.plot(tarr_odd * 24., oddmod/1.e-6, color='#e88b1c', lw=4)
        self.ax4.axhline(y=odd_depth/1.e-6, color='b', ls='--', lw=2)
        self.ax4.get_yaxis().get_major_formatter().set_useOffset(False)
        ymin = np.min(np.min(f_odd), np.min(f_even))
        self.ax4.set_ylim(ymin/1.e-6, 6*sigma/1.e-6)

        self.ax5.plot(dt_even * 24., f_even/1.e-6, lw=0, marker='.', color='0.7')
        self.ax5.plot(dt_even_bin * 24., f_even_bin/1e-6, marker='.', color='#1c79e8', markersize=20, lw=0)
        self.ax5.plot(tarr_even * 24., evenmod/1.e-6, color='#e88b1c', lw=4)
        self.ax5.axhline(y=even_depth/1.e-6, color='b', ls='--', lw=2)
        self.ax5.annotate(r'Diff: %.3f $\sigma$' % difference, xy=(0.52, 0.1), xycoords='figure fraction', size=15)
        self.ax5.get_yaxis().get_major_formatter().set_useOffset(False)
        self.ax5.set_ylim(ymin/1.e-6, 6*sigma/1.e-6)
        self.ax5.get_yaxis().get_major_formatter().set_useOffset(False)

   
    
    
    def plot_figure6to7(self):
        self.ax6.set_xlabel("Time from T0 (hours)")
        self.ax6.set_ylabel("Delta X (Pixel)")
    
        self.ax6.annotate("Centroids", [0.71, 0.1], xycoords='figure fraction',  fontweight='black', color= '#be7013', fontsize=18, fontstyle='oblique')
    
       
        self.ax7.set_ylabel('Delta Y(Pixel)', color='r')
        self.ax7.tick_params('y', colors='r')
        indexlist = self.transit.get_intran(self.jd)
        x_clean = np.array([])
        y_clean = np.array([])
        dt_clean = np.array([])
        import scipy as sp
        from scipy import interpolate
        for i in xrange(len(indexlist)):

            #print x
            x_now = self.x[indexlist[i]]
            y_now = self.y[indexlist[i]]
            dt_now = self.transit.dt[indexlist[i]] 

            if len(x_now) < 5:
                fit = np.polyfit(dt_now, x_now, 1)
                xmod = np.polyval(fit, dt_now)

                fit = np. polyfit(dt_now, y_now, 1)
                ymod = np.polyval(fit, dt_now)

            else:
                tck = sp.interpolate.splrep(dt_now, x_now, s=len(dt_now) - np.sqrt(2 * len(dt_now)))
                xmod = sp.interpolate.splev(dt_now, tck)
                sigx = np.std(x_now - xmod)

                tck = sp.interpolate.splrep(dt_now, y_now, s=len(dt_now) - np.sqrt(2 * len(dt_now)))
                ymod = sp.interpolate.splev(dt_now, tck)
                sigy = np.std(y_now - ymod)

                try:
                    # sigma clip
                    good_x = np.where(abs(x_now-xmod) < 2*sigx)[0]
                    good_y = np.where(abs(y_now-ymod) < 2*sigy)[0]

                    tck = sp.interpolate.splrep(dt_now[good_x], x_now[good_x], s=len(good_x) - np.sqrt(2 * len(good_x)))
                    xmod = sp.interpolate.splev(dt_now, tck)

                    tck = sp.interpolate.splrep(dt_now[good_y], y_now[good_y], s=len(good_y) - np.sqrt(2 * len(good_y)))
                    ymod = sp.interpolate.splev(dt_now, tck)

                except TypeError:
                    pass
            # plt.plot(t_now, x_now, 'k.')
            # plt.plot(t_now, xmod, color='r', lw=2)
            # plt.show()

            x_now -= xmod
            y_now -= ymod

            x_clean = np.append(x_clean, x_now)
            y_clean = np.append(y_clean, y_now)
            dt_clean = np.append(dt_clean, dt_now)

        sig = np.std(x_clean)
        y_clean += max(x_clean) - min(y_clean) + 0.8*sig
        self.ax6.plot(dt_clean*24., x_clean, 'k.')
        self.ax7.plot(dt_clean*24., y_clean, 'r.')

        self.ax6.set_ylim(min(x_clean), max(y_clean))  # there's probably an easier way to set the same ylim for the two axes
        self.ax7.set_ylim(min(x_clean), max(y_clean))
 
        return


