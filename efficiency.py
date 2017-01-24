import glob
import math
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from datetime import datetime
from joblib import Parallel, delayed
import multiprocessing


matplotlib.use('TkAgg')
startTime = datetime.now()

count = 0
counter = 0
meansum = []
#splines = []
#fluxes = []
#start = 196009 					#KIC 12557548 = 196009
#finish = 196010
#start = 81083 						#KIC 7303287 = 81083
#finish = 81084
start = 0
finish = 10
allfiles = glob.glob('/media/marko/Saul/Detrenddata1/*.dat')


#for filename in allfiles[start:finish]:
def plotting(filename):
	#if count%20 == 0: print count
	destination = "Efficiency/" + filename[31:43]
	tempFile = open(filename, 'r')
	temptime = []
	tempflux = []
	time = []
	flux = []

	headerInfo = tempFile.readline()
	radius = headerInfo.split()[2]
	temperature = headerInfo.split()[3]
	magnitude = headerInfo.split()[4]
	namestr = filename[31:43]
	radstr = "Radius = " + radius
	tempstr = "Temp = " + temperature
	magstr = "Mag = " + magnitude

	tempFile.next()
	tempFile.next()

	for line in tempFile:
		parts = line.split()
		temptime.append(float(parts[0]))
		tempflux.append(float(parts[1]))

	time = temptime[0:len(time)-1]
	flux = tempflux[0:len(flux)-1]
	meansum.append(sum(flux)/len(flux))

	ti = time[0]
	tf = time[len(time)-1]
	tdelta = time[1]-time[0]
	nbins = (tf-ti)/tdelta
	for n in range(0, 20):
		if pow(2,n) > nbins:
			nbins = pow(2,n)
			break
	fftarray = np.zeros(int(nbins)*8, dtype=float)
	boxcarfftarray = [1] * nbins*8

	#print ti, tf, tdelta, fftarray, len(time)

	for value in range(0, len(time)-1):
		binint = (time[value] - time[0])/(time[1]-time[0])
#		if flux[value] > 5*rms:
#			continue
		fftarray[int(binint)] = flux[value]-1
		boxcarfftarray[int(binint)] = flux[value]

#	smoothedarray = np.zeros(int(nbins)*8, dtype=float)
	smoothedarray = [1] * nbins*8
#	for mainbin in range(0, 29):
	for mainbin in range(0, len(fftarray)-1):
		boxsum = []
		binsused = 0
		if mainbin-24 < 0 or mainbin+24 > len(fftarray)-1:
			binmin = 0
			binmax = mainbin
		else:
			binmin = mainbin-24
			binmax = mainbin+24
		for sidebins in range(binmin, binmax):
			if boxcarfftarray[sidebins] != 1:
				boxsum.append(boxcarfftarray[sidebins])
				binsused += 1
			else:
				continue
		if binsused == 0:
			continue
		smoothedarray[mainbin] = ((sum(boxsum)/binsused))

#	subbedarray = np.subtract(boxcarfftarray, smoothedarray)
	subbedarray = np.zeros(len(boxcarfftarray), dtype=float)
	for x in range(0, len(subbedarray)):
		if boxcarfftarray[x] == 1:
			subbedarray[x] = 0
		else:
			subbedarray[x] = boxcarfftarray[x] - smoothedarray[x]
	binsarray = np.linspace(0, len(subbedarray)-1, len(subbedarray))
	binsarray = binsarray*tdelta

#	print fftarray[0:99], smoothedarray[0:99], subbedarray[0:99]
	#(time - time0)/tdelta = indicator of which integer it falls under = which bin it falls under
	fftarrayout = np.fft.fft(fftarray)
	fftarrayabs = []
	for bin in range(0, len(fftarrayout)-1):
		fftarrayabs.append(abs(fftarrayout[bin]))

	fftlen = np.linspace(0, len(fftarrayabs), len(fftarrayabs))
	freqarray = (fftlen-1)/((tdelta)*8*nbins)

#-----------------FFT of Smoothed Boxcar Subtracted Data---------------------------
	fftarrayout2 = np.fft.fft(subbedarray)
	fftarrayabs2 = []
	for bin in range(0, len(fftarrayout2)-1):
		fftarrayabs2.append(abs(fftarrayout2[bin]))

	fftlen2 = np.linspace(0, len(fftarrayabs2), len(fftarrayabs2))
	freqarray2 = (fftlen2-1)/((tdelta)*8*nbins)
#-----------------FFT of Smoothed Boxcar Subtracted Data---------------------------

#	fig = plt.figure()
#	sp1 = fig.add_subplot(111)
#	sp1.plot(time, flux, 'k')
#	plt.figtext(0.755, 0.865, namestr, fontsize = 'small')
#	plt.figtext(0.755, 0.845, radstr[0:len(radstr)-3], fontsize = 'small')
#	plt.figtext(0.755, 0.825, tempstr[0:len(tempstr)-7], fontsize = 'small')
#	plt.figtext(0.755, 0.805, magstr[0:len(magstr)-3], fontsize = 'small')
#	sp1.get_yaxis().get_major_formatter().set_useOffset(False)
#	plt.savefig(destination)
#	plt.close()
#	count += 1
	return

num_cores = multiprocessing.cpu_count()
results = Parallel(n_jobs=num_cores)(delayed(plotting)(filename) for filename in allfiles[start:finish])
endTime = datetime.now()
print "Time for program to run = ", endTime-startTime