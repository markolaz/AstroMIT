#112 samples per orbit
#start at JD 2073, line 17857
#end at JD 2076.09, line 17857 + (112*40) = 22337
#do the same for the second outburst

import matplotlib.pyplot as plt
import numpy as np

def modulus(x, y):
	return float(x)/float(y) - int(x/y)

rawtime = []
rawflux = []

with open('GaussElimData/V893-ap2.dat') as infile:
	for line in infile:
		rawtime.append(float(line.split()[0]))
		rawflux.append(float(line.split()[1]))

period = 0.07596192

time1, flux1 = rawtime[17857:22337], rawflux[17857:22337]

time2, flux2 = rawtime[64559:70159], rawflux[64559:70159]

times1 = [time1[i*112:(i+1)*112] for i in xrange(int(len(time1)/112))]
fluxes1 = [flux1[i*112:(i+1)*112] for i in xrange(int(len(flux1)/112))]

times2 = [time2[i*112:(i+1)*112] for i in xrange(int(len(time2)/112))]
fluxes2 = [flux2[i*112:(i+1)*112] for i in xrange(int(len(flux2)/112))]

orbitwindows1 = [0] * len(times1)
orbitwindows2 = [0] * len(times2)

count = 0

for timewindow, fluxwindow in zip(times1, fluxes1):
	foldarr = [0] * 100
	exp = [0] * 100
	for timevalue, fluxvalue in zip(timewindow, fluxwindow):
		fraction = modulus(float(timevalue), period)
		foldarr[int(fraction*100)] += float(fluxvalue)/np.mean(fluxwindow)
		exp[int(fraction*100)] += 1
	folded = np.divide(foldarr, exp)
	orbitwindows1[count] = folded + float(count)/10
	count += 1
	del foldarr, exp, fraction, folded

count = 0

for timewindow, fluxwindow in zip(times2, fluxes2):
	foldarr = [0] * 100
	exp = [0] * 100
	for timevalue, fluxvalue in zip(timewindow, fluxwindow):
		fraction = modulus(float(timevalue), period)
		foldarr[int(fraction*100)] += float(fluxvalue)/np.mean(fluxwindow)
		exp[int(fraction*100)] += 1
	folded = np.divide(foldarr, exp)
	orbitwindows2[count] = folded + float(count)/10
	count += 1
	del foldarr, exp, fraction, folded

phase = np.linspace(-0.5, 0.5, 100)

plt.figure(figsize=(8,16))
plt.title('First Flare Orbit Fold')
for curve in range(len(orbitwindows1)):
	plt.step(phase, orbitwindows1[curve], 'k')
ax = plt.gca()
ax.set_xlim([-0.51, 0.51])
plt.xlabel('Orbital Phase')
plt.ylabel('Flux, incremented by 0.1 for each orbit')
plt.savefig('V893 Sco/FirstFlareOrbits.png')

plt.close()

plt.figure(figsize=(8,16))
plt.title('Second Flare Orbit Fold')
for curve in range(len(orbitwindows2)):
	plt.step(phase, orbitwindows2[curve], 'k')
ax = plt.gca()
ax.set_xlim([-0.51, 0.51])
plt.xlabel('Orbital Phase')
plt.ylabel('Flux, incremented by 0.1 for each orbit')
plt.savefig('V893 Sco/SecondFlareOrbits.png')