import matplotlib.pyplot as plt
import numpy as np
from math import *
from datetime import datetime

def FFT(datafile, time=[], flux=[], oversampling=8, nbins=None, sample_rate=None, plotting=False):

	from collections import Counter

	if not (time and flux):
		time, flux = [], []
		with open(datafile) as infile:
			for line in infile:
				time.append(float(line.split()[0]))
				flux.append(float(line.split()[1]))

		infile.close()

	if not sample_rate: # Sample rate has units of time (sec), not frequency (Hz)
		sample_rates = [a - b for a, b in zip(time, time[1:])] 
		mode = Counter(sample_rates)
		sample_rate = abs(mode.most_common(1)[0][0]) # Gets the mode, or most common, sample rate and uses that for the FFT
		print sample_rate

	if not nbins:
		for n in range(0, 30):
			if pow(2, n) > ((time[-1]-time[0])/sample_rate):
				nbins = pow(2, n+1)
	
	fft_input = np.zeros(int(nbins)*oversampling, dtype=float)

	fluxmean = np.mean(flux)
	
	for value in range(len(time)):
		binint = (float(time[value]) - float(time[0]))/(sample_rate)
		fft_input[int(binint)] = flux[value]-fluxmean
	
	fft_output = np.fft.fft(fft_input)

	fft = map(abs, fft_output)
	
	fftlen = np.linspace(0, len(fft), len(fft))
	freq = (fftlen-1)/(sample_rate*oversampling*nbins)

	if plotting:
		fig = plt.figure()
		plt.plot(freq, fft, 'k')
		ax = plt.gca()
		ax.set_xlim([0, 1./(2*sample_rate)])
		return fig, freq, fft
	else:
		return freq, fft

data = FFT('GaussElimData/V893-ap2_filtered_div.dat', oversampling=8, nbins=pow(2, 18), plotting=False)
plt.show()
