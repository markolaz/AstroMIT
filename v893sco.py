import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from math import *

def modulus(x, y):
	return float(x)/float(y) - int(x/y)
def mean(array):
	return sum(array)/len(array)

rawtime = []
rawflux = []
time = []
flux = []

count = 0

with open('GaussElimData/V893-ap2.dat') as infile:
	for line in infile:
		if line.split()[0] < 2170. or (line.split()[0] > 2180. and line.split()[0] < 2200.) or line.split()[0] > 2210.:
			rawtime.append(float(line.split()[0]))
			rawflux.append(float(line.split()[1]))

with open('GaussElimData/V893-ap2_filtered_div.dat') as infile:
	for line in infile:
		if line.split()[0] < 2170. or (line.split()[0] > 2180. and line.split()[0] < 2200.) or line.split()[0] > 2210.:
			time.append(float(line.split()[0]))
			flux.append(float(line.split()[1]))

period = 0.07596192

#-----------------Full Data Set Fold-------------------------

foldarray = [0] * 200
exposure = [0] * 200
for i in range(len(time)):
	frac = modulus((float(time[i])+period/3.5), period)
	foldarray[int(frac*200)] += float(flux[i])
	exposure[int(frac*200)] += 1
fold = np.divide(foldarray, exposure)
phase = np.linspace(-0.5, 0.5, 200)

days = [time[i*1470:(i+1)*1470] for i in xrange(int(len(time)/1470))]
fluxes = [flux[i*1470:(i+1)*1470] for i in xrange(int(len(flux)/1470))]

daywindows = [0] * len(days)

#-----------------1 Day Folds----------------------------

for timewindow, fluxwindow in zip(days, fluxes):
	foldarr = [0] * 200
	exp = [0] * 200
	for timevalue, fluxvalue in zip(timewindow, fluxwindow):
		fraction = modulus((float(timevalue)+period/3.5), period)
		foldarr[int(fraction*200)] += float(fluxvalue)
		exp[int(fraction*200)] += 1
	folded = np.divide(foldarr, exp)
	daywindows[count] = folded
	count += 1
	del foldarr, exp, fraction, folded

#-----------------Chi2 Goodness of Fit with scaling parameters alpha & beta to fit full dataset fold to 1 day folds-------------------------

alpha = []

diffgrid = np.array([[0.]*231 for _ in range(200)])
diffgrid.shape = (231, 200)

for value in range(len(daywindows)):

	fmatrix = np.array([[0.]*2 for _ in range(2)])
	fmatrix.shape = (2, 2)
	dmatrix = np.array([0.]*2)
	dmatrix.shape = (2,1)
	
	fmatrix[0, 0] = sum([x**2 for x in fold])
	fmatrix[1, 0] = sum([x*1 for x in fold])
	fmatrix[0, 1] = sum([x*1 for x in fold])
	fmatrix[1, 1] = sum([1 for x in fold])
	
	dmatrix[0, 0] = sum([x*y for x, y in zip(fold, daywindows[value])])
	dmatrix[1, 0] = sum([y*1 for y in daywindows[value]])
	
	params = np.dot(np.linalg.inv(fmatrix), dmatrix)

	chi2 = daywindows[value] - fold*params[0] - len(daywindows[value])*params[1]

	time_str = '{}{}'.format('Time = ', time[value*1470])

	diff = daywindows[value] - fold

	for bin in range(len(fold)):
		diffgrid[int(value*3), bin] = diff[bin]
		diffgrid[int(value*3+1), bin] = diff[bin]
		diffgrid[int(value*3+2), bin] = diff[bin]

	ax = plt.gca()
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	ax.get_yaxis().get_major_formatter().set_useOffset(False)
	plt.step(phase, (fold*params[0] + params[1]), 'red')
	plt.step(phase, daywindows[value], 'k')
	plt.figtext(0.655, 0.165, time_str, fontsize = 'medium')
	plt.xlabel('Orbital Phase')
	plt.ylabel('Normalized Flux')
	ax.set_ylim([0.3, 1.3])
	figname = "{}%02d.png".format('V893 Sco/chisquareplot') % value
	plt.savefig(figname)
	plt.close()

	alpha.append(float(params[0]))

alphabins = np.linspace(0, len(alpha)-1, len(alpha))

#-----------------Dividing the original filtered (by way of division and gaussian elimination) dataset by the fold based on fractional relationship to the fold-----------------

flux_div = np.zeros_like(flux)
for value in range(len(flux_div)):
	fraction = modulus((float(time[value])+period/3.5), period)
	bin_flux = fold[int(fraction*200)]
	flux_div[value] = float(flux[value])/bin_flux

nptime = np.array(time, dtype='float')
npflux = np.array(flux_div, dtype='float')
data = np.array([nptime, npflux])
data = data.T
with open('V893 Sco/V893_Sco_divfilt_divfold.txt', 'w') as infile:
	np.savetxt(infile, data, fmt=['%f', '%f'])

#-----------------FFT of raw data-----------------

nbins = pow(2, 18)
tdelta = float(time[1]) - float(time[0])

fftarray = np.zeros(int(nbins)*8, dtype=float)

rawfluxmean = mean(rawflux)

for value in range(len(rawtime)):
	binint = (float(rawtime[value]) - float(rawtime[0]))/(float(rawtime[1])-float(rawtime[0]))
	fftarray[int(binint)] = rawflux[value]-rawfluxmean

fftarrayout = np.fft.fft(fftarray)
fftarrayabs = []
for bin in range(0, len(fftarrayout)-1):
	fftarrayabs.append(abs(fftarrayout[bin]))

fftlen = np.linspace(0, len(fftarrayabs), len(fftarrayabs))
freqarray = (fftlen-1)/((tdelta)*8*nbins)

#-----------------FFT of divided original filtered data-----------------

fftarrayraw = np.zeros(int(nbins)*8, dtype=float)
fftarrayfilt = np.zeros(int(nbins)*8, dtype=float)

fluxmean = mean(flux)
div_fluxmean = mean(flux_div)

for value in range(len(time)):
	binint = (float(time[value]) - float(time[0]))/(float(time[1])-float(time[0]))
	fftarrayraw[int(binint)] = flux[value]-fluxmean
	fftarrayfilt[int(binint)] = flux_div[value]-div_fluxmean

fftarrayoutraw = np.fft.fft(fftarrayraw)
fftarrayoutfilt = np.fft.fft(fftarrayfilt)
fftarrayabsraw = []
fftarrayabsfilt = []
for bin in range(0, len(fftarrayoutraw)-1):
	fftarrayabsraw.append(abs(fftarrayoutraw[bin]))
	fftarrayabsfilt.append(abs(fftarrayoutfilt[bin]))

#-----------------Folding over period acquired from FFT-----------------
offset = 7000
peakfreq = np.argmax(fftarrayabsraw[offset:(len(fftarrayabsraw)/2)])
period_div = 1/(float(freqarray[peakfreq+offset]))*2
foldarray_div = [0] * 200
exposure_div = [0] * 200
for i in range(len(flux_div)):
	frac_div = modulus((float(time[i])+period_div/5), period_div)
	foldarray_div[int(frac_div*200)] += float(flux_div[i])
	exposure_div[int(frac_div*200)] += 1
fold_div = np.divide(foldarray_div, exposure_div)
phase_div = np.linspace(-0.5, 0.5, 200)
period_str = '{}{}'.format('Period = ', period_div)

#--------------Plots------------------

#--------------Plot of fold of division filtered data divided by overall data fold--------------

plt.step(phase_div, fold_div, 'k')
plt.savefig('V893 Sco/filtereddividedfold.png')
plt.close()

#--------------Plot of division filtered data divided by the overall data fold--------------

plt.plot(time, flux_div, 'k')
ax1 = plt.gca()
ax1.set_xlim([int(float(time[0])), int(float(time[-1]))])
plt.savefig('V893 Sco/dividedbyfold.png')
plt.close()

#--------------Plot of FFT of raw data--------------

plt.plot(freqarray, fftarrayabs, 'k')
ax2 = plt.gca()
ax2.set_xlim([0, 734])
ax2.set_ylim([0, 7e7])
plt.figtext(0.600, 0.755, period_str, fontsize='medium')
plt.savefig('V893 Sco/raw_FFT.png')
#plt.show()
plt.close()

#--------------Plot of FFT of division filtered data--------------

plt.plot(freqarray, fftarrayabsraw, 'k')
ax2 = plt.gca()
ax2.set_xlim([0, 734])
ax2.set_ylim([0, 4000])
plt.figtext(0.600, 0.755, period_str, fontsize='medium')
plt.savefig('V893 Sco/filt_FFT.png')
#plt.show()
plt.close()

#--------------Plot of FFT of division filtered data divided by overall data fold--------------

plt.plot(freqarray, fftarrayabsfilt, 'k')
ax2 = plt.gca()
ax2.set_xlim([0, 734])
ax2.set_ylim([0, 4000])
plt.figtext(0.600, 0.755, period_str, fontsize='medium')
plt.savefig('V893 Sco/filt_divfold_FFT.png')
plt.close()

#--------------Plot of scale factor vs day (77 days)--------------

plt.step(alphabins, alpha, 'k')
plt.ylabel('Scale Factor Value')
plt.xlabel('Day')
plt.savefig('V893 Sco/scalefactorevolution.png')
plt.close()

#--------------Plot of differences between 1 day folds and the template for all 77 days with 3 bins per day--------------

plt.imshow(diffgrid, cmap='hot')
plt.colorbar()
plt.xlabel('Bins')
plt.ylabel('Days (in groups of 3)')
plt.title('Difference between Template and 1 Day Folds')
plt.savefig('V893 Sco/Template-1Day-diff.png')
hdu = fits.PrimaryHDU(diffgrid)
hdu.writeto('V893 Sco/Template-1Day-diff.fits', overwrite=True)
plt.close()




#sum(data)*T - alpha*sum(T)*T - beta*n*T = 0

#sum(data)*T = alpha*sum(T)*T + beta*n*T

#go back through fold, take the flux from the fold bin found by finding the phase, divide data by fold, make new data series