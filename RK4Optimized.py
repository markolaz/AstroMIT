import numpy as np
from datetime import datetime

def vectors(n_particles, v_ratio):

	anglesrand = np.random.rand(n_particles, 2)

	thetavec = np.arccos(1. - 2.*anglesrand[:, 0], dtype='float32')
	phivec = 2.*np.pi*anglesrand[:, 1]
	
	uvector = v_ratio*np.sin(thetavec)*np.cos(phivec)+1./np.sqrt(2)
	vvector = v_ratio*np.sin(thetavec)*np.sin(phivec)+1./np.sqrt(2)
	wvector = v_ratio*np.cos(thetavec)

	return thetavec, phivec, uvector, vvector, wvector

def RK4(dt, beta, x, y, z, u, v, w):

    ix_0 = dt*(u)                                                         #dt*f(t, x, u)
    jx_0 = dt*(-x*(1-beta)/(np.sqrt(x**2+y**2+z**2))**3)                            #dt*g(t, x, u)

    iy_0 = dt*(v)                                                         #dt*f(t, x, u)
    jy_0 = dt*(-y*(1-beta)/(np.sqrt(x**2+y**2+z**2))**3)                                                #dt*g(t, x, u)

    iz_0 = dt*(w)                           
    jz_0 = dt*(-z*(1-beta)/(np.sqrt(x**2+y**2+z**2))**3)
    
    ix_1 = dt*(u + jx_0/2.)                                               #dt*f(t + dt/2., x + ix_0/2., u + jx_0/2.)
    jx_1 = dt*(-(x + ix_0/2.)*(1.-beta)/(np.sqrt((x + ix_0/2.)**2+(y + iy_0/2.)**2+(z + iz_0/2.)**2))**3)                        #dt*g(t + dt/2., x + ix_0/2., u + jx_0/2.)
    
    iy_1 = dt*(v + jy_0/2.)                                               #dt*f(t + dt/2., x + ix_0/2., u + jx_0/2.)
    jy_1 = dt*(-(y + iy_0/2.)*(1.-beta)/(np.sqrt((x + ix_0/2.)**2+(y + iy_0/2.)**2+(z + iz_0/2.)**2))**3)                        #dt*g(t + dt/2., x + ix_0/2., u + jx_0/2.)

    iz_1 = dt*(w + jz_0/2.)                                               
    jz_1 = dt*(-(z + iz_0/2.)*(1.-beta)/(np.sqrt((x + ix_0/2.)**2+(y + iy_0/2.)**2+(z + iz_0/2.)**2))**3)
    
    ix_2 = dt*(u + jx_1/2.)                                               #dt*f(t + dt/2., x + ix_1/2., u + jx_1/2.)
    jx_2 = dt*(-(x + ix_1/2.)*(1.-beta)/(np.sqrt((x + ix_1/2.)**2+(y + iy_1/2.)**2+(z + iz_1/2.)**2))**3)                        #dt*g(t + dt/2., x + ix_1/2., u + jx_1/2.)
    
    iy_2 = dt*(v + jy_1/2.)                                               #dt*f(t + dt/2., x + ix_1/2., u + jx_1/2.)
    jy_2 = dt*(-(y + iy_1/2.)*(1.-beta)/(np.sqrt((x + ix_1/2.)**2+(y + iy_1/2.)**2+(z + iz_1/2.)**2))**3)                        #dt*g(t + dt/2., x + ix_1/2., u + jx_1/2.)

    iz_2 = dt*(w + jz_1/2.)                                               #dt*f(
    jz_2 = dt*(-(z + iz_1/2.)*(1.-beta)/(np.sqrt((x + ix_1/2.)**2+(y + iy_1/2.)**2+(z + iz_1/2.)**2))**3)
    
    ix_3 = dt*(u + jx_2)                                                  #dt*f(t + dt, x + ix_2, u + jx_2)
    jx_3 = dt*(-(x + ix_2)*(1.-beta)/(np.sqrt((x + ix_2)**2+(y + iy_2)**2+(z + iz_2)**2))**3)                              #dt*g(t + dt, x + ix_2, u + jx_2)
    
    iy_3 = dt*(v + jy_2)                                                  #dt*f(t + dt, x + ix_2, u + jx_2)
    jy_3 = dt*(-(y + iy_2)*(1.-beta)/(np.sqrt((x + ix_2)**2+(y + iy_2)**2+(z + iz_2)**2))**3)                              #dt*g(t + dt, x + ix_2, u + jx_2)

    iz_3 = dt*(w + jz_2)                                               
    jz_3 = dt*(-(z + iz_2)*(1.-beta)/(np.sqrt((x + ix_2)**2+(y + iy_2)**2+(z + iz_2)**2))**3)
    
    x = x + (1./6.)*(ix_0 + 2*ix_1 + 2*ix_2 + ix_3)
    u = u + (1./6.)*(jx_0 + 2*jx_1 + 2*jx_2 + jx_3)
    
    y = y + (1./6.)*(iy_0 + 2*iy_1 + 2*iy_2 + iy_3)
    v = v + (1./6.)*(jy_0 + 2*jy_1 + 2*jy_2 + jy_3)

    z = z + (1./6.)*(iz_0 + 2*iz_1 + 2*iz_2 + iz_3)
    w = w + (1./6.)*(jz_0 + 2*jz_1 + 2*jz_2 + jz_3)

    return x, y, z, u, v, w

#x'' = -GMx/sqrt(x^2+y^2)^3
#y'' = -GMy/sqrt(x^2+y^2)^3
#
#x' = u                     1st eq, f(t, x, y)
#u' = -GMx/(x^2+y^2)^3      2nd eq, g(t, x, y)
#
#y' = v
#v' = -GMy/(x^2+y^2)^3

def main():

	start = datetime.now()

	n_particles = 1000  # 4.5 minutes for 100 particles, also 4.5 minutes for 1000 particles
	
	t = 0.0
	
	x0 = 1./np.sqrt(2)
	y0 = -1./np.sqrt(2)
	z0 = 0.0
	
	v_ratio = 0.008 # the v/V ratio, where v is v_particle and V is v_body
	beta0 = 0.
	
	dt = 0.1
	
	# impact parameter, b, is the z-coordinate
	# D is distance from particle (center of ribbon) to center of star
	# L is full distance of ribbon projection
	
	xvector = np.full(n_particles, x0, 'float32')
	yvector = np.full(n_particles, y0, 'float32')
	zvector = np.full(n_particles, z0, 'float32')
	
	uvector = np.zeros(n_particles, 'float32')
	vvector = np.zeros(n_particles, 'float32')
	wvector = np.zeros(n_particles, 'float32')
	
	xrotvec = np.zeros(n_particles, 'float32')
	yrotvec = np.zeros(n_particles, 'float32')
	zrotvec = np.zeros(n_particles, 'float32')
	
	betavec = np.full(n_particles, beta0, 'float32')
	
	thetavec, phivec, uvector, vvector, wvector = vectors(n_particles, v_ratio)
	
	timespan = np.linspace(0, 199.9, 2000) #end time = (dt*n_steps)-dt, or (0.1*1000)-0.1
	
	beta = betavec
	
	for t in timespan:
		xvector, yvector, zvector, uvector, vvector, wvector = RK4(dt, beta, xvector, yvector, zvector, uvector, vvector, wvector)
	
	print datetime.now()-start

if __name__ == '__main__':
	main()
