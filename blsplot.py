import matplotlib.pyplot as plt
import numpy as np

#KIC 8644288

#infile = open('/home/marko/Desktop/MIT Research/FFITools/LCTools/etc/data/raw/12351927.dat')
#lines = infile.readlines()
#lines = [str(line.split()[1] + '\t' + str(2.0-float(line.split()[2])) + '\t1\n') for line in lines]

#with open('/home/marko/Desktop/MIT Research/FFITools/LCTools/etc/data/raw/12351927_err.dat', 'w') as outfile:
#	outfile.writelines(lines)
#	outfile.close()

period = []
sn = []
time = []
flux = []

with open('/home/marko/Desktop/MIT Research/FFITools/LCTools/etc/data/raw/12351927_err.dat.bls') as infile:
	infile.next()
	for line in infile:
		period.append(float(line.split()[0]))
		sn.append(float(line.split()[1]))


with open('/home/marko/Desktop/MIT Research/FFITools/LCTools/etc/data/raw/12351927.dat') as infile:
	infile.next()
	for line in infile:
		time.append(float(line.split()[1]))
		flux.append(float(line.split()[2]))

plt.plot(period, sn, 'k')
plt.xlabel('Period (Days)')
plt.ylabel('Signal-to-Noise')
plt.show()