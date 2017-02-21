import matplotlib.pyplot as plt
import numpy as np
import math

n_params = 21

#file = open('/home/marko/Desktop/MIT Research/KIC012557548.dat', 'r')
file = open('/home/marko/Desktop/MIT Research/KIC008153499.dat', 'r')
file.next()
file.next()

time = []
flux = []

data_period = 25.    #in days
shortest_period = 2.5 #in days
f_max = 1./shortest_period
delta_f = 1./data_period

freqrange = np.linspace(delta_f, f_max, f_max/delta_f)

for line in file:
	parts = line.split()
	time.append(float(parts[0]))
	flux.append(float(parts[1]))

maxtime = int(data_period*round(float(time[-2])/data_period))

timesplit = np.linspace(data_period, maxtime, maxtime/data_period)

#timefull = np.zeros_like(time)
#fitfull = np.zeros_like(flux)
#datafull = np.zeros_like(flux)

timefull = []
fitfull = []
datafull = []

for val in xrange(len(timesplit)):
#for val in range(0, len(timesplit)):
	if val==maxtime: continue
	time50 = [t for t in time if timesplit[val] <= t <= timesplit[val+1]]
	data50 = [flux[time.index(t)] for t in time50]

	fmatrix = np.array([[0.]*n_params for _ in range(n_params)])
	fmatrix.shape = (n_params,n_params)
	dmatrix = np.array([0.]*n_params)
	dmatrix.shape = (n_params,1)

	if time50==[]:continue
	
	for x in xrange(len(time50)-1):
	#	d = (19.2 + 16.5*(t/100.) - 134*((t/100.)**2) + 449*((t/100.)**3))
		t = time50[x]
		d = data50[x]
		f1 = math.sin(2*math.pi*freqrange[0]*t)
		f2 = math.sin(2*math.pi*freqrange[1]*t)
		f3 = math.sin(2*math.pi*freqrange[2]*t)
		f4 = math.sin(2*math.pi*freqrange[3]*t)
		f5 = math.sin(2*math.pi*freqrange[4]*t)
		f6 = math.sin(2*math.pi*freqrange[5]*t)
		f7 = math.sin(2*math.pi*freqrange[6]*t)
		f8 = math.sin(2*math.pi*freqrange[7]*t)
		f9 = math.sin(2*math.pi*freqrange[8]*t)
		f10 = math.sin(2*math.pi*freqrange[9]*t)
		f11 = math.cos(2*math.pi*freqrange[0]*t)
		f12 = math.cos(2*math.pi*freqrange[1]*t)
		f13 = math.cos(2*math.pi*freqrange[2]*t)
		f14 = math.cos(2*math.pi*freqrange[3]*t)
		f15 = math.cos(2*math.pi*freqrange[4]*t)
		f16 = math.cos(2*math.pi*freqrange[5]*t)
		f17 = math.cos(2*math.pi*freqrange[6]*t)
		f18 = math.cos(2*math.pi*freqrange[7]*t)
		f19 = math.cos(2*math.pi*freqrange[8]*t)
		f20 = math.cos(2*math.pi*freqrange[9]*t)
		f21 = 1.
	
		fvector = [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, f19, f20, f21]
	
		for row in range(n_params):
			for column in range(n_params):
				fmatrix[row, column] += (fvector[row]*fvector[column])
			dmatrix[row] += d*fvector[row]
	
	
	invfmatrix = np.linalg.inv(fmatrix)
	
	solmatrix = np.dot(invfmatrix, dmatrix)
	
	#print solmatrix
	
	fit50 = []
	
	for x in xrange(len(time50)):
		t = time50[x]
		fit1 =  solmatrix[0]*math.sin(2*math.pi*freqrange[0]*t)
		fit2 =  solmatrix[1]*math.sin(2*math.pi*freqrange[1]*t)
		fit3 =  solmatrix[2]*math.sin(2*math.pi*freqrange[2]*t)
		fit4 =  solmatrix[3]*math.sin(2*math.pi*freqrange[3]*t)
		fit5 =  solmatrix[4]*math.sin(2*math.pi*freqrange[4]*t)
		fit6 =  solmatrix[5]*math.sin(2*math.pi*freqrange[5]*t)
		fit7 =  solmatrix[6]*math.sin(2*math.pi*freqrange[6]*t)
		fit8 =  solmatrix[7]*math.sin(2*math.pi*freqrange[7]*t)
		fit9 =  solmatrix[8]*math.sin(2*math.pi*freqrange[8]*t)
		fit10 = solmatrix[9]*math.sin(2*math.pi*freqrange[9]*t)
		fit11 = solmatrix[10]*math.cos(2*math.pi*freqrange[0]*t)
		fit12 = solmatrix[11]*math.cos(2*math.pi*freqrange[1]*t)
		fit13 = solmatrix[12]*math.cos(2*math.pi*freqrange[2]*t)
		fit14 = solmatrix[13]*math.cos(2*math.pi*freqrange[3]*t)
		fit15 = solmatrix[14]*math.cos(2*math.pi*freqrange[4]*t)
		fit16 = solmatrix[15]*math.cos(2*math.pi*freqrange[5]*t)
		fit17 = solmatrix[16]*math.cos(2*math.pi*freqrange[6]*t)
		fit18 = solmatrix[17]*math.cos(2*math.pi*freqrange[7]*t)
		fit19 = solmatrix[18]*math.cos(2*math.pi*freqrange[8]*t)
		fit20 = solmatrix[19]*math.cos(2*math.pi*freqrange[9]*t)
		fit21 = solmatrix[20]*1.
	
		fit50.extend(data50[x]-(fit1+fit2+fit3+fit4+fit5+fit6+fit7+fit8+fit9+fit10+fit11+fit12+fit13+fit14+fit15+fit16+fit17+fit18+fit19+fit20+fit21))

	#timefull.extend(time50)
	datafull.extend(data50)
	fitfull.extend(fit50)
	#for x in xrange(len(time50)):
	#	fitfull[((val-1)*len(time50))+x] = fit50[x]
	#	timefull[((val-1)*len(time50))+x] = time50[x]
	#	datafull[((val-1)*len(time50))+x] = data50[x]
	#fig = plt.figure()
	#sp1 = fig.add_subplot(211)
	#sp1.step(time50, data50, 'k')
	#
	#sp2 = fig.add_subplot(212)
	#sp2.step(time50, fit50, 'k')
#
	#figname = 'fig' + str(val)
	#plt.savefig(figname)
#
	#plt.close()

	del time50, fit50, fvector, data50, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, f19, f20, fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9, fit10, fit11, fit12, fit13, fit14, fit15, fit16, fit17, fit18, fit19, fit20

#fig2 = plt.figure()
#
#splt1 = fig2.add_subplot(211)
#splt1.step(time[:-2], flux[:-2], 'k')
#
#splt2 = fig2.add_subplot(212)
#splt2.step(timefull, fitfull, 'k')
#
#plt.show()











#timevector = []
#datavector = []

#f11vector = []
#f12vector = []
#f13vector = []
#f14vector = []
#f21vector = []
#f22vector = []
#f23vector = []
#f24vector = []
#f31vector = []
#f32vector = []
#f33vector = []
#f34vector = []
#f41vector = []
#f42vector = []
#f43vector = []
#f44vector = []

#d1vector = []
#d2vector = []
#d3vector = []
#d4vector = []

#for t in xrange(100):
#	d = (19.2 + 16.5*(t/100.) - 134*((t/100.)**2) + 449*((t/100.)**3))	
	#f11vector.append(1*1)
	#f12vector.append(1*(t/100.))
	#f13vector.append(1*((t/100.)**2))
	#f14vector.append(1*((t/100.)**3))
	#f21vector.append(1*(t/100.))
	#f22vector.append((t/100.)*(t/100.))
	#f23vector.append((t/100.)*((t/100.)**2))
	#f24vector.append((t/100.)*((t/100.)**3))
	#f31vector.append(1*((t/100.)**2))
	#f32vector.append((t/100.)*((t/100.)**2))
	#f33vector.append(((t/100.)**2)*((t/100.)**2))
	#f34vector.append(((t/100.)**2)*((t/100.)**3))
	#f41vector.append(1*((t/100.)**3))
	#f42vector.append((t/100.)*((t/100.)**3))
	#f43vector.append(((t/100.)**2)*((t/100.)**3))
	#f44vector.append(((t/100.)**3)*((t/100.)**3))
	#d1vector.append(d*1)
	#d2vectoar.append(d*(t/100.))
	#d3vector.append(d*((t/100.)**2))
	#d4vector.append(d*((t/100.)**3))


##	timevector.append(t)
##	datavector.append(d)

#f11 = sum(f11vector)
#f12 = sum(f12vector)
#f13 = sum(f13vector)
#f14 = sum(f14vector)
#f21 = sum(f21vector)
#f22 = sum(f22vector)
#f23 = sum(f23vector)
#f24 = sum(f24vector)
#f31 = sum(f31vector)
#f32 = sum(f32vector)
#f33 = sum(f33vector)
#f34 = sum(f34vector)
#f41 = sum(f41vector)
#f42 = sum(f42vector)
#f43 = sum(f43vector)
#f44 = sum(f44vector)
#d1 = sum(d1vector)
#d2 = sum(d2vector)
#d3 = sum(d3vector)
#d4 = sum(d4vector)
#

#dmatrix = np.matrix([[d1], [d2], [d3], [d4]])
#fmatrix = np.matrix([[f11, f12, f13, f14], [f21, f22, f23, f24], [f31, f32, f33, f34], [f41, f42, f43, f44]])
#invfmatrix = np.linalg.inv(fmatrix)

#solmatrix = np.dot(invfmatrix, dmatrix)

#print solmatrix