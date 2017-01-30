import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.interpolate as inter
import glob
from joblib import Parallel, delayed
import multiprocessing
from datetime import datetime

#!/usr/bin/env python -W ignore::VisibleDeprecationWarning

startTime = datetime.now()

def remove_sv(desig,spline_width):
    print '-----------------------'+desig+'-----------------------'
    if desig == 'KIC008153499': return
    if desig == 'KIC008153500': return
    if desig == 'KIC008153533': return
    if desig == 'KIC008153547': return
    '''
    file = 'data/'+desig+'_raw.txt'
    newdata = pd.read_table(file, sep='\s+', header = 0, names=["Time", "Flux", "unc","quarter"])
    

    #print newdata
    time = np.array(newdata["Time"], 'd')
    flux = np.array(newdata["Flux"], 'd')
    unc = np.array(newdata["unc"], 'd')
    quarter = np.array(newdata["quarter"], 'd')
    '''
#    file = 'data/'+desig+'_flux.txt'
#    if len(desig) == 8: file = 'data/raw/KIC0'+desig+'.dat'
#    if len(desig) == 7: file = 'data/raw/KIC00'+desig+'.dat'
    file = '/media/marko/Saul/Detrenddata1/'+desig+'.dat'
    newdata = pd.read_table(file, sep='\s+', header = 2, names=["Time", "Flux", "Quarter"])
    time_all = np.array(newdata["Time"], 'd')
    flux_all = np.array(newdata["Flux"], 'd')
    Quarter = np.array(newdata["Quarter"], 'i')
    unc = np.ones_like(flux_all)*np.std(flux_all[0:20])
    #remove the transits first from further analysis
    
    
    #define the effective quarter number i.e. when there is a data gap, increase the quarter number
    #Quarter = np.zeros_like(time_all)
    
    qua = 0
    n_gap =len(time_all)-1
    for i in range(n_gap):
        dif = time_all[i+1]-time_all[i]
        Quarter[i] = qua
        if dif > 30./60./24.0*2:
            #print dif
            qua += 1
    Quarter[-1] = qua

#print np.unique(Quarter)

    text_file = open("data/splined/"+desig+"_detrended.txt", "w")
#    text_file = open("data/"+desig+"_detrended_5days.txt", "w")
    for i in range(np.max(Quarter)):
        #print 'quarter', i+1
        time = time_all[Quarter == i]
        flux = flux_all[Quarter == i]
        flux_notran = flux#[(abs((time-t0) % p) >= 0.75*dura) & (abs((time-t0) % p) <= (p-0.75*dura))]
    
        unc_notran = unc#[(abs((time-t0) % p) >= 0.75*dura) & (abs((time-t0) % p) <= (p-0.75*dura))]
        time_notran = time#[(abs((time-t0) % p) >= 0.75*dura) & (abs((time-t0) % p) <= (p-0.75*dura))]
        '''
        plt.scatter(time % p,flux, color = 'green')
        plt.scatter(time_notran % p,flux_notran)
        plt.show()
        
        plt.scatter(time,flux, color = 'green')
        plt.scatter(time_notran,flux_notran)
        plt.show()
        '''
    #detrend light curve using splines of 0.75 days duration, iterate a few times
        if len(time_notran) <= 1: continue
        time_fit = time_notran
        flux_fit = flux_notran
        flux_spline = flux_fit
    
    #find the gaps in time
        time_dif = time_fit[1:-1]-time_fit[0:-2]
        time_dif2 = np.array([0.]+list(time_dif)) # the other side of gap
        gaps = list(time_fit[time_dif > 30./60./24.0*2])+list(time_fit[time_dif2 > 30./60./24.0*2])
        
        
        #plt.scatter(time_notran,flux_notran)
        #plt.scatter(time_clean,flux_clean, color = 'green')
        #plt.scatter(gaps,np.ones_like(gaps),color = 'red')
        #plt.plot(time_notran,ss(time_notran),color = 'red')
        #plt.show()
        
        #spline_width = p
    
        breakpoint = np.linspace(np.min(time_notran),np.max(time_notran),int((np.max(time_notran)-np.min(time_notran))/spline_width))
        breakpoint = list(breakpoint)+list(gaps)
        #print np.min(time_fit),np.max(time_fit),np.min(breakpoint),np.max(breakpoint)
        #print gaps
        breakpoint = np.sort(np.array(list(breakpoint)))
    
        for m in range(5):
            diff_std = np.std(flux_fit-flux_spline)
            time_fit = time_fit[abs(flux_fit-flux_spline)<=(5.*diff_std)]
            flux_fit = flux_fit[abs(flux_fit-flux_spline)<=(5.*diff_std)]

        #print breakpoint[1:-1],np.min(time_fit),np.max(time_fit)
        #print breakpoint
            ss = inter.LSQUnivariateSpline(time_fit,flux_fit,breakpoint[1:-1])
            flux_spline = ss(time_fit)
#print len([abs(flux_fit-flux_spline)<=(5.*diff_std)])

#produce the spline detrended light curve
        flux_spl=flux_notran/ss(time_notran)

        flux_spline = ss(time_notran)

#remove the outlier
        diff_std = np.std(flux_notran-flux_spline)
        time_clean = time_notran[abs(flux_notran-flux_spline)<=(5.*diff_std)]
        flux_clean = flux_notran[abs(flux_notran-flux_spline)<=(5.*diff_std)]


        '''
        plt.scatter(time_notran,flux_notran)
        plt.scatter(time_clean,flux_clean, color = 'green')
        plt.scatter(gaps,np.ones_like(gaps),color = 'red')
        plt.plot(time_notran,ss(time_notran),color = 'red')
        #plt.show()
        '''
    
    #let's see if the spline successful remove stellar variability
        flux_spline = ss(time_clean)
#plt.scatter((time_clean - t0) ,flux_clean/flux_spline)
        #plt.show()
    
    #finally plot the phase folded light curve and bin them in phase
    
#plt.scatter(((time_clean - t0) % p)/p ,flux_clean/flux_spline)
#time_bin,flux_bin,unc_bin = bin((time_clean - t0) % p, flux_clean/flux_spline, nbins,p)
#plt.errorbar(time_bin/p,flux_bin,yerr= unc_bin, color = 'orange')
        #plt.show()
        #print "total", len(time_notran), "clean", len(time_clean)


#now output the phase-binned light curve for further analysis
#text_file = open("data/"+desig+"_phase_bin.txt", "w")
        for j in range(len(time_clean)):
            text_file.write(str(time_clean[j])+'   ')
            text_file.write(str(flux_clean[j]/flux_spline[j])+'\n')
                #text_file.write(str(unc_bin[i])+'\n')
    text_file.close()

def bin(time, flux, nbins,p):
    time_mod = np.ones_like(time)
    t_min = np.min(time)
    length = p
    time_mod[:] = time[:]#-t_min
    time_bin = []
    flux_bin = []
    unc_bin = []
    for i in range(nbins):
        inside = [abs(time_mod - length*(i+0.5)/(nbins*1.0)) < length/(nbins*2.0)]
        
        if len(flux[inside])>1:
            time_bin.append(np.mean(time_mod[inside]))#+t_min
            flux_bin.append(np.mean(flux[inside]))
            tmpp= np.sqrt(len(flux[inside])*1.)
            unc_bin.append(np.std(flux)/tmpp)
        elif len(flux[inside]) == 1:
            time_bin.append(np.mean(time_mod[inside]))
            flux_bin.append(np.mean(flux[inside]))
            unc_bin.append(np.nanstd(flux))
    
    return np.array(time_bin),np.array(flux_bin),np.array(unc_bin)

rawdatafiles = glob.glob('/media/marko/Saul/Detrenddata1/*.dat')


start = 100500      
finish = 101000   #if starting at 0, make finish the exact number you want, not n-1


splinewidth = 5.0        #default 1.5


#num_cores = multiprocessing.cpu_count()
results = Parallel(n_jobs=-1)(delayed(remove_sv)(filename[31:43], splinewidth) for filename in rawdatafiles[start:finish])
#for filename in rawdatafiles[start:finish]:
#    if filename[31:43] == 'KIC008153499': continue
#    results = Parallel(n_jobs=-1)(delayed(remove_sv)(filename[31:43], splinewidth))
print 'Time to run = ', datetime.now() - startTime

#remove_sv('72.01',sap =0, llc = 1)#note sap ==1 use sap flux, sap ==0 use PDC flux; if llc ==0 use short cadence; if slc ==1 use long cadence