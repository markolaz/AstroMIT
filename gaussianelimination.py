import matplotlib.pyplot as plt
import numpy as np
import math

n_params = 4

fmatrix = np.array([[0.]*n_params for _ in range(n_params)])
fmatrix.shape = (n_params,n_params)
dmatrix = np.array([0.]*n_params)
dmatrix.shape = (n_params,1)

for t in xrange(100):
	d = (19.2 + 16.5*(t/100.) - 134*((t/100.)**2) + 449*((t/100.)**3))

	f1 = 1
	f2 = (t/100.)
	f3 = (t/100.)**2
	f4 = (t/100.)**3

	fvector = [f1, f2, f3, f4]

	for row in range(n_params):
		for column in range(n_params):
			fmatrix[row, column] += (fvector[row]*fvector[column])
		dmatrix[row] += d*fvector[row]

invfmatrix = np.linalg.inv(fmatrix)

solmatrix = np.dot(invfmatrix, dmatrix)

print solmatrix














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
	#d2vector.append(d*(t/100.))
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