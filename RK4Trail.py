import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm
from astropy.io import fits
import math
import numpy as np
import datetime

beginning_time = datetime.datetime.now()

#x'' = -GMx/sqrt(x^2+y^2+z^2)^3
#y'' = -GMy/sqrt(x^2+y^2+z^2)^3
#z'' = -GMz/sqrt(x^2+y^2+z^2)^3
#
#x' = u                         1st eq, f(t, x, y)
#u' = -GMx/(x^2+y^2+z^2)^3      2nd eq, g(t, x, y)
#
#y' = v
#v' = -GMy/(x^2+y^2+z^2)^3
#
#z' = w
#w' = -GMz/(x^2+y^2+z^2)^3


n_particles = 10000  #1 minute for 100 particles, linear scaling (10x more particles takes 10x longer, only for trail)
n_steps = 25000
dt = 0.02
tau = 2616.5 #decay constant for particles to decay once they are 180 degrees out of phase with the parent body

start_time = 0
end_time = (dt*n_steps)-dt


t = 0.0

x0 = 0.0
y0 = -1.
z0 = 0.0

xvector = [x0]*n_particles
yvector = [y0]*n_particles
zvector = [z0]*n_particles

uvector = [0]*n_particles
vvector = [0]*n_particles
wvector = [0]*n_particles

xrotvec = [0]*n_particles
yrotvec = [0]*n_particles
zrotvec = [0]*n_particles

v_particle = 0.01 #actually is the v/V ratio, where v is v_particle and V is v_body
beta = 0.

anglesrand = np.random.rand(n_particles, 2)
betaran = np.random.rand(n_particles)

for particle in xrange(n_particles):

    theta = math.acos(1. - 2.*anglesrand[particle, 0])
    phi = 2.*math.pi*anglesrand[particle, 1]

    u0 = v_particle*math.sin(theta)*math.cos(phi)+1
    v0 = v_particle*math.sin(theta)*math.sin(phi)
    w0 = v_particle*math.cos(theta)

    uvector[particle] = u0
    vvector[particle] = v0
    wvector[particle] = w0

#v_ratiostr = "{} = {}".format("v/V", v_particle)
#betastr = "{} = {}".format("Beta", beta)

#real units of time are acquired by using T * P / 2pi, where T is t in loop and P is period of 1 radian in orbit

timespan = np.linspace(start_time, end_time, n_steps) #end time = (dt*n_steps)-dt, or (0.1*1000)-0.1

xarray = np.array([[0.]*n_particles for _ in range(n_steps)])
yarray = np.array([[0.]*n_particles for _ in range(n_steps)])
zarray = np.array([[0.]*n_particles for _ in range(n_steps)])

WDgrid = np.array([[0.]*601 for _ in range(601)])
WDgrid.shape = (601, 601)

for t in timespan:

    if t%100 == 0: print t, datetime.datetime.now()-beginning_time

    for particle in xrange(n_particles):

        row = int(t*10)
        column = particle

        x1 = xvector[particle]
        y1 = yvector[particle]
        z1 = zvector[particle]

        u = uvector[particle]
        v = vvector[particle]
        w = wvector[particle]

        ix_0 = dt*(u)                                                         #dt*f(t, x1, u)
        jx_0 = dt*(-x1*(1-beta)/(math.sqrt(x1**2+y1**2+z1**2))**3)                            #dt*g(t, x1, u)
    
        iy_0 = dt*(v)                                                         #dt*f(t, x1, u)
        jy_0 = dt*(-y1*(1-beta)/(math.sqrt(x1**2+y1**2+z1**2))**3)                                                #dt*g(t, x1, u)

        iz_0 = dt*(w)                           
        jz_0 = dt*(-z1*(1-beta)/(math.sqrt(x1**2+y1**2+z1**2))**3)
        
        ix_1 = dt*(u + jx_0/2.)                                               #dt*f(t + dt/2., x1 + ix_0/2., u + jx_0/2.)
        jx_1 = dt*(-(x1 + ix_0/2.)*(1.-beta)/(math.sqrt((x1 + ix_0/2.)**2+(y1 + iy_0/2.)**2+(z1 + iz_0/2.)**2))**3)                        #dt*g(t + dt/2., x1 + ix_0/2., u + jx_0/2.)
    
        iy_1 = dt*(v + jy_0/2.)                                               #dt*f(t + dt/2., x1 + ix_0/2., u + jx_0/2.)
        jy_1 = dt*(-(y1 + iy_0/2.)*(1.-beta)/(math.sqrt((x1 + ix_0/2.)**2+(y1 + iy_0/2.)**2+(z1 + iz_0/2.)**2))**3)                        #dt*g(t + dt/2., x1 + ix_0/2., u + jx_0/2.)

        iz_1 = dt*(w + jz_0/2.)                                               
        jz_1 = dt*(-(z1 + iz_0/2.)*(1.-beta)/(math.sqrt((x1 + ix_0/2.)**2+(y1 + iy_0/2.)**2+(z1 + iz_0/2.)**2))**3)
        
        ix_2 = dt*(u + jx_1/2.)                                               #dt*f(t + dt/2., x1 + ix_1/2., u + jx_1/2.)
        jx_2 = dt*(-(x1 + ix_1/2.)*(1.-beta)/(math.sqrt((x1 + ix_1/2.)**2+(y1 + iy_1/2.)**2+(z1 + iz_1/2.)**2))**3)                        #dt*g(t + dt/2., x1 + ix_1/2., u + jx_1/2.)
    
        iy_2 = dt*(v + jy_1/2.)                                               #dt*f(t + dt/2., x1 + ix_1/2., u + jx_1/2.)
        jy_2 = dt*(-(y1 + iy_1/2.)*(1.-beta)/(math.sqrt((x1 + ix_1/2.)**2+(y1 + iy_1/2.)**2+(z1 + iz_1/2.)**2))**3)                        #dt*g(t + dt/2., x1 + ix_1/2., u + jx_1/2.)

        iz_2 = dt*(w + jz_1/2.)                                               #dt*f(
        jz_2 = dt*(-(z1 + iz_1/2.)*(1.-beta)/(math.sqrt((x1 + ix_1/2.)**2+(y1 + iy_1/2.)**2+(z1 + iz_1/2.)**2))**3)
        
        ix_3 = dt*(u + jx_2)                                                  #dt*f(t + dt, x1 + ix_2, u + jx_2)
        jx_3 = dt*(-(x1 + ix_2)*(1.-beta)/(math.sqrt((x1 + ix_2)**2+(y1 + iy_2)**2+(z1 + iz_2)**2))**3)                              #dt*g(t + dt, x1 + ix_2, u + jx_2)
    
        iy_3 = dt*(v + jy_2)                                                  #dt*f(t + dt, x1 + ix_2, u + jx_2)
        jy_3 = dt*(-(y1 + iy_2)*(1.-beta)/(math.sqrt((x1 + ix_2)**2+(y1 + iy_2)**2+(z1 + iz_2)**2))**3)                              #dt*g(t + dt, x1 + ix_2, u + jx_2)

        iz_3 = dt*(w + jz_2)                                               
        jz_3 = dt*(-(z1 + iz_2)*(1.-beta)/(math.sqrt((x1 + ix_2)**2+(y1 + iy_2)**2+(z1 + iz_2)**2))**3)
    
        x2 = x1 + (1./6.)*(ix_0 + 2*ix_1 + 2*ix_2 + ix_3)
        u = u + (1./6.)*(jx_0 + 2*jx_1 + 2*jx_2 + jx_3)
    
        y2 = y1 + (1./6.)*(iy_0 + 2*iy_1 + 2*iy_2 + iy_3)
        v = v + (1./6.)*(jy_0 + 2*jy_1 + 2*jy_2 + jy_3)

        z2 = z1 + (1./6.)*(iz_0 + 2*iz_1 + 2*iz_2 + iz_3)
        w = w + (1./6.)*(jz_0 + 2*jz_1 + 2*jz_2 + jz_3)

        sine = math.sin(t)
        cosine = math.cos(t)

        xvector[particle] = x2
        yvector[particle] = y2
        zvector[particle] = z2

        uvector[particle] = u
        vvector[particle] = v
        wvector[particle] = w

        xrotvec[particle] = x2*cosine + y2*sine
        yrotvec[particle] = y2*cosine - x2*sine
        zrotvec[particle] = z2

        xarray[row, column] = x1*cosine + y1*sine
        yarray[row, column] = y1*cosine - x1*sine
        zarray[row, column] = z1

        if abs(xarray[row, column]) <= 0.03 and abs(zarray[row, column]) <= 0.03:

            x_grid = int(round(xarray[row, column]*10000+300, 0))
            z_grid = int(round(zarray[row, column]*10000+300, 0))

            WDgrid[z_grid, x_grid] += 5*math.exp(-t/tau)

    t = t + dt

hdu = fits.PrimaryHDU(WDgrid)
hdu.writeto('ParticleDensity.fits', overwrite=True)

print 'Time taken = ', datetime.datetime.now()-beginning_time
