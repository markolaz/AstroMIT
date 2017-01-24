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
#start = 196009 					#KIC 12557548 = 196009
#finish = 196010
#start = 81083 						#KIC 7303287 = 81083
#finish = 81084
start = 0
finish = 1000
rawdatafiles = glob.glob('/media/marko/Saul/Detrenddata1/*.dat')
#splinedfiles = glob.glob('/home/marko/Desktop/MIT Research/SplineFit/data/splined/*.txt')


def BoxCarSubFFT(filename):
	#print (filename[31:43])
	destination = "Efficiency/Raw&SplineFFT/" + filename[31:43]
#	if len(filename) == 69: destination = "Efficiency/" + filename[48:55]
#	if len(filename) == 70: destination = "Efficiency/" + filename[48:56]
#	destination = "KIC Output Images/Interesting/" + filename[31:43] + "logpart"
#	destination = "KIC Output Images/Interesting/" + filename[31:43] + "lin"


	rawFile = open(filename, 'r')
	splFile = open('/home/marko/Desktop/MIT Research/SplineFit/data/splined/'+filename[31:43]+'_detrended.txt', 'r')
	temptimeraw = []
	tempfluxraw = []
	temptimespl = []
	tempfluxspl = []
	rawTime = []
	rawFlux = []
	splTime = []
	splFlux = []


	headerInfo = rawFile.readline()
	radius = headerInfo.split()[2]
	temperature = headerInfo.split()[3]
	magnitude = headerInfo.split()[4]
	namestr = filename[31:43]
	radstr = "Radius = " + radius
	tempstr = "Temp = " + temperature
	magstr = "Mag = " + magnitude

	rawFile.next()
	rawFile.next()

	for line in rawFile:
		parts = line.split()
		temptimeraw.append(float(parts[0]))
		tempfluxraw.append(float(parts[1]))
	for line in splFile:
		temptimespl.append(float(line.split()[0]))
		tempfluxspl.append(float(line.split()[1]))


	rawTime = temptimeraw[0:len(rawTime)-1]
	rawFlux = tempfluxraw[0:len(rawFlux)-1]
	splTime = temptimespl[0:len(splTime)-1]
	splFlux = tempfluxspl[0:len(splFlux)-1]
	meansum.append(sum(rawFlux)/len(rawFlux))

#	rms = np.std(rawFlux)

	ti = rawTime[0]
	tf = rawTime[len(rawTime)-1]
	tdelta = rawTime[1]-rawTime[0]
	nbins = (tf-ti)/tdelta 	#(rawTime - time0)/tdelta = indicator of which integer it falls under = which bin it falls under
	for n in range(0, 20):
		if pow(2,n) > nbins:
			nbins = pow(2,n)
			break
	fftarrayraw = np.zeros(int(nbins)*8, dtype=float)
	fftarrayspl = np.zeros(int(nbins)*8, dtype=float)

	#print ti, tf, tdelta, fftarrayraw, len(rawTime)

	for value in range(0, len(rawTime)-1):
		binint = (rawTime[value] - rawTime[0])/(rawTime[1]-rawTime[0])
		fftarrayraw[int(binint)] = rawFlux[value]-1
	for value in range(0, len(splTime)-1):
		binint = (splTime[value] - splTime[0])/(splTime[1]-splTime[0])
		fftarrayspl[int(binint)] = splFlux[value]-1

	binsarrayraw = np.linspace(0, len(fftarrayraw)-1, len(fftarrayraw))
	binsarrayraw = binsarrayraw*tdelta
	binsarrayspl = np.linspace(0, len(fftarrayspl)-1, len(fftarrayspl))
	binsarrayspl = binsarrayspl*tdelta

	fftarrayoutraw = np.fft.fft(fftarrayraw)
	fftarrayoutspl = np.fft.fft(fftarrayspl)
	fftarrayabsraw = []
	fftarrayabsspl = []
	for bin in range(0, len(fftarrayoutraw)-1):
		fftarrayabsraw.append(abs(fftarrayoutraw[bin]))
	for bin in range(0, len(fftarrayoutspl)-1):
		fftarrayabsspl.append(abs(fftarrayoutspl[bin]))
#	[(fftarrayabsraw.append(abs(fftarrayoutraw[bin])) for bin in range(len(fftarrayoutraw)))]
#	mean = np.mean(fftarrayabsraw)
	fftlenraw = np.linspace(0, len(fftarrayabsraw), len(fftarrayabsraw))
	fftlenspl = np.linspace(0, len(fftarrayabsspl), len(fftarrayabsspl))
	freqarrayraw = (fftlenraw-1)/((tdelta)*8*nbins)
	freqarrayspl = (fftlenspl-1)/((tdelta)*8*nbins)



#	output = open('Efficiency/'+filename[31:43]+'FFT.txt', 'w+')
#	data = np.array([freqarrayraw[0:len(freqarrayraw)/2], fftarrayabsraw[0:len(fftarrayabsraw)/2]])
#	data = data.T
#	np.savetxt(output, data, fmt=['%.6f','%.6f'])



	fig = plt.figure()
	fig.set_size_inches(12, 8)
	sp1 = fig.add_subplot(221)
	sp1.plot(rawTime, rawFlux, 'k')
	plt.figtext(0.755/2, 0.865, namestr, fontsize = 'small')
	plt.figtext(0.755/2, 0.845, radstr[0:len(radstr)-3], fontsize = 'small')
	plt.figtext(0.755/2, 0.825, tempstr[0:len(tempstr)-7], fontsize = 'small')
	plt.figtext(0.755/2, 0.805, magstr[0:len(magstr)-3], fontsize = 'small')
	sp1.get_yaxis().get_major_formatter().set_useOffset(False) #eliminates exponential representation of values (1.06, not 9.990^-1 + 0.061)
	sp1.set_title("Raw Data")
	#sp1.set_ylim([0.95, 1.05])
	sp2 = fig.add_subplot(222)
	sp2.plot(splTime, splFlux, 'k')
	plt.figtext(0.755+0.045, 0.865, namestr, fontsize = 'small')
	plt.figtext(0.755+0.045, 0.845, radstr[0:len(radstr)-3], fontsize = 'small')
	plt.figtext(0.755+0.045, 0.825, tempstr[0:len(tempstr)-7], fontsize = 'small')
	plt.figtext(0.755+0.045, 0.805, magstr[0:len(magstr)-3], fontsize = 'small')
	sp2.get_yaxis().get_major_formatter().set_useOffset(False) #eliminates exponential representation of values (1.06, not 9.990^-1 + 0.061)
	sp2.set_title("Filtered Data")
	sp3 = fig.add_subplot(223)
	sp3.plot(freqarrayraw[0:len(freqarrayraw)/2], fftarrayabsraw[0:len(fftarrayabsraw)/2], 'k')
	plt.figtext(0.755/2, 0.865-0.865/2, namestr, fontsize = 'small')
	plt.figtext(0.755/2, 0.845-0.865/2, radstr[0:len(radstr)-3], fontsize = 'small')
	plt.figtext(0.755/2, 0.825-0.865/2, tempstr[0:len(tempstr)-7], fontsize = 'small')
	plt.figtext(0.755/2, 0.805-0.865/2, magstr[0:len(magstr)-3], fontsize = 'small')
	sp3.get_yaxis().get_major_formatter().set_useOffset(False) #eliminates exponential representation of values (1.06, not 9.990^-1 + 0.061)
	sp3.set_xlim([0, 25])
	sp3.set_title("Raw Data FFT")
	sp4 = fig.add_subplot(224)
	sp4.plot(freqarrayspl[0:len(freqarrayspl)/2], fftarrayabsspl[0:len(fftarrayabsspl)/2], 'k')
	plt.figtext(0.755+0.045, 0.865-0.865/2, namestr, fontsize = 'small')
	plt.figtext(0.755+0.045, 0.845-0.865/2, radstr[0:len(radstr)-3], fontsize = 'small')
	plt.figtext(0.755+0.045, 0.825-0.865/2, tempstr[0:len(tempstr)-7], fontsize = 'small')
	plt.figtext(0.755+0.045, 0.805-0.865/2, magstr[0:len(magstr)-3], fontsize = 'small')
	sp4.get_yaxis().get_major_formatter().set_useOffset(False) #eliminates exponential representation of values (1.06, not 9.990^-1 + 0.061)
	sp4.set_xlim([0, 25])
	sp4.set_title("Filtered Data FFT")
	plt.savefig(destination)
	plt.close()
#	sp2.set_yscale('log')
#	sp2.set_ylim([0.3, 100])		#for log part, no limits for log full
#	sp2.set_ylim([0, 30])			# for linear

	del rawFile, splFile, temptimeraw, temptimespl, tempfluxraw, tempfluxspl, rawTime, splTime, rawFlux, splFlux, parts#, fig
	return

num_cores = multiprocessing.cpu_count()
results = Parallel(n_jobs=num_cores)(delayed(BoxCarSubFFT)(filename) for filename in rawdatafiles[start:finish])

print "rawTime for program to run = ", datetime.now() - startTime