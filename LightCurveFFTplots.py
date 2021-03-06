import glob
import math
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from datetime import datetime


matplotlib.use('TkAgg')
startTime = datetime.now()

count = 0
#counter = 0
meansum = []
#splines = []
#fluxes = []
start = 196009 					#KIC12557548 = 196009
finish = 196010
#start = 0
#finish = 10
allfiles = glob.glob('/media/marko/Saul/Detrenddata1/*.dat')

#for filename in allfiles:				#specific star finder
#	if filename == '/media/marko/Saul/Detrenddata1/KIC012557548.dat':
#		print counter
#	counter+=1


for filename in allfiles[start:finish]:
	#print (filename[31:43])
	print count
#	destination = "KIC Output Images/" + filename[31:43]
	destination = "KIC Output Images/Interesting/" + filename[31:43] + "logpart"


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
#	rms = np.std(flux)

	ti = time[0]
	tf = time[len(time)-1]
	tdelta = time[1]-time[0]
	nbins = (tf-ti)/tdelta
	for n in range(0, 20):
		if pow(2,n) > nbins:
			nbins = pow(2,n)
			break
	fftarray = np.zeros(int(nbins)*8, dtype=float)

	#print ti, tf, tdelta, fftarray, len(time)

	for value in range(0, len(time)-1):
		binint = (time[value] - time[0])/(time[1]-time[0])
#		if flux[value] > 5*rms:
#			continue
		fftarray[int(binint)] = flux[value]-1
	#(time - time0)/tdelta = indicator of which integer it falls under = which bin it falls under
	fftarrayout = np.fft.fft(fftarray)
	fftarrayabs = []
	for bin in range(0, len(fftarrayout)-1):
		fftarrayabs.append(abs(fftarrayout[bin]))

	fftlen = np.linspace(0, len(fftarrayabs), len(fftarrayabs))
	freqarray = (fftlen-1)/((tdelta)*8*nbins)

	fig = plt.figure()
	sp1 = fig.add_subplot(211)
	sp1.plot(time, flux, 'k')
	plt.figtext(0.755, 0.865, namestr, fontsize = 'small')
	plt.figtext(0.755, 0.845, radstr[0:len(radstr)-3], fontsize = 'small')
	plt.figtext(0.755, 0.825, tempstr[0:len(tempstr)-7], fontsize = 'small')
	plt.figtext(0.755, 0.805, magstr[0:len(magstr)-3], fontsize = 'small')
	sp1.get_yaxis().get_major_formatter().set_useOffset(False) #eliminates exponential representation of values (1.06, not 9.990^-1 + 0.061)
	#sp1.set_ylim([0.95, 1.05])
	sp2 = fig.add_subplot(212)
	sp2.plot(freqarray[0:len(freqarray)/2], fftarrayabs[0:len(fftarrayabs)/2], 'k')
	sp2.set_xlim([0, 25])
	sp2.set_yscale('log')
	sp2.set_ylim([0.3, 100])		#for log part, no limits for log full
#	sp2.set_ylim([0, 30])		# for linear
	plt.savefig(destination)
	plt.close()

	count+=1
	if count==(finish-start)/2:
		midTime = datetime.now()
	del tempFile, temptime, tempflux, time, flux, parts

avg = sum(meansum)/count
print("Average Flux per Star = %f" % avg)

print "First Half Time = ", (midTime - startTime)
print "Second Half Time = ", (datetime.now() - midTime)



#spline = interpolate.splrep(time, flux, s=(len(time)-math.sqrt(2*len(time))))
#spline = interpolate.splrep(time, flux, s=(20))
#fluxnew = interpolate.splev(time, spline)
#splines.append(spline)
#fluxes.append(fluxnew)
#subbed = fluxnew - flux


#plt.figure()
#plt.hist(meansum, 10000, [0.98, 1.02])
#plt.draw()
#plt.savefig("KICAvgFluxDistribution")
#plt.close()