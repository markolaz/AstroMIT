import matplotlib.pyplot as plt
import numpy as np
import math

beta = 0.0005
v_ratio = 0.005
n_particles = 300
anglesrand = np.random.rand(n_particles, 2)
a = 1.
p = 4000.

for t in xrange(1000):

	xvector = []
	yvector = []
	zvector = []
	xrotvec = []
	yrotvec = []
	zrotvec = []

	for particle in xrange(n_particles):

		theta = math.acos(1. - 2.*anglesrand[particle, 0])
		phi = 2.*math.pi*anglesrand[particle, 1]
	
		xi = 2.*v_ratio*math.sin(theta)*math.cos(phi) + v_ratio**2
		chi = 2.*v_ratio*math.sin(theta)*math.cos(phi) + (v_ratio**2)*(math.sin(theta)**2*math.cos(phi)**2 + math.cos(phi)**2)
	
		new_a = a*(1.-beta)/(1.-2.*beta-xi)
		mu = math.atan(v_ratio*math.cos(theta)/(1. + v_ratio*math.sin(theta)*math.cos(phi)))
	
		new_p = p*math.pow(new_a/a, 3./2.)/math.sqrt(1.-beta)
	
	
		x0 = a*math.sin(2.*math.pi*(t*40.)/p)
		y0 = -a*math.cos(2.*math.pi*(t*40.)/p)
		z0 = 0.

		X = new_a*math.sin(2*math.pi*(t*40.)/new_p)*math.cos(mu)
		Y = -new_a*math.cos(2*math.pi*(t*40.)/new_p)
		Z = new_a*math.sin(2*math.pi*(t*40.)/new_p)*math.sin(mu)

		phi2 = math.atan2(X, -Y) #azimuthal angle used for removing the rotation
		phi0 = math.atan2(x0, -y0)
		delta_phi = phi2-phi0

		x_rotframe = new_a*math.sin(delta_phi)
		y_rotframe = -new_a*math.cos(delta_phi)
		z_rotframe = Z

		x_frame = X - x0
		y_frame = Y - y0
		z_frame = Z - z0

		if not ((x_frame**2+z_frame**2) < (0.011)**2 and y_frame > 0):	
			xvector.append(x_frame)
			yvector.append(y_frame)
			zvector.append(z_frame)

		xrotvec.append(x_rotframe)
		yrotvec.append(y_rotframe)
		zrotvec.append(z_rotframe)


	

	xzphase = "{} = {}".format("Phase", (t*40)/4000.)
	v_ratiostr = "{} = {}".format("v/V", v_ratio)
	betastr = "{} = {}".format("Beta", beta)

	circle1 = plt.Circle((0,0), 0.011, color='cyan')
	circle2 = plt.Circle((0,0), 0.011, color='cyan')
	circle3 = plt.Circle((0,0), 0.011, color='cyan')

	fig1 = plt.figure()
	xzaxes = plt.gca()
	xzaxes.set_xlim([-0.6, 0.6])
	xzaxes.set_ylim([-0.5, 0.5])
	xzaxes.add_artist(circle1)
	if t%100. == 0: plt.plot(xvector, zvector, 'go', ms=1)
	else: plt.plot(xvector, zvector, 'ro', ms=1)
	plt.figtext(0.705, 0.865, xzphase, fontsize = 'large')
	plt.figtext(0.705, 0.825, v_ratiostr, fontsize = 'large')
	plt.figtext(0.705, 0.785, betastr, fontsize = 'large')
	xzfigname = '{}{}'.format('Exploding Orbit Images/xz_image', t)
	fig1.suptitle('XZ Plane W/ Rot')
	plt.savefig(xzfigname)
	plt.close()

	fig2 = plt.figure()
	xyaxes = plt.gca()
	xyaxes.set_xlim([-1.2, 1.2])
	xyaxes.set_ylim([-1.2, 1.2])
	xyaxes.add_artist(circle2)
	plt.plot(xrotvec, yrotvec, 'bo', ms=1)
	plt.figtext(0.705, 0.825, v_ratiostr, fontsize = 'large')
	plt.figtext(0.705, 0.785, betastr, fontsize = 'large')
	xyfigname = '{}{}'.format('Exploding Orbit Images/xy_image', t)
	fig2.suptitle('XY Plane')
	plt.savefig(xyfigname)
	plt.close()

	fig3 = plt.figure()
	xzaxesrot = plt.gca()
	xzaxesrot.set_xlim([-0.6, 0.6])
	xzaxesrot.set_ylim([-0.5, 0.5])
	xzaxesrot.add_artist(circle3)
	if t%100. == 0: plt.plot(xrotvec, zrotvec, 'go', ms=1)
	else: plt.plot(xrotvec, zrotvec, 'ro', ms=1)
	plt.figtext(0.705, 0.865, xzphase, fontsize = 'large')
	plt.figtext(0.705, 0.825, v_ratiostr, fontsize = 'large')
	plt.figtext(0.705, 0.785, betastr, fontsize = 'large')
	xzfignamerot = '{}{}'.format('Exploding Orbit Images/xz_imagerot', t)
	fig1.suptitle('XZ Plane No Rot')
	plt.savefig(xzfignamerot)
	plt.close()

#	if t%100. == 0:
#		xzaxes2 = plt.gca()
#		xzaxes2.set_xlim([-0.8, 0.8])
#		xzaxes2.set_ylim([-0.01, 0.01])
#		plt.plot(xvector, zvector, 'ro')
#		xzfigname2 = '{}{}'.format('Exploding Orbit Images/xzfront_image', t)
#		fig1.suptitle('XZ Plane')
#		plt.savefig(xzfigname2)
#		plt.close()	

	del fig1, fig2, fig3



#print anglesrand, anglesrand[2,0], anglesrand[2,1]