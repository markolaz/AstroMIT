import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np
from datetime import datetime

#x'' = -GMx/sqrt(x^2+y^2)^3
#y'' = -GMy/sqrt(x^2+y^2)^3
#
#x' = u                     1st eq, f(t, x, y)
#u' = -GMx/(x^2+y^2)^3      2nd eq, g(t, x, y)
#
#y' = v
#v' = -GMy/(x^2+y^2)^3

start = datetime.now()

n_particles = 1000  #4.5 minutes for 100 particles, also 4.5 minutes for 1000 particles

t = 0.0

x0 = 1./math.sqrt(2)
y0 = -1./math.sqrt(2)
z0 = 0.0

# impact parameter, b, is the z-coordinate
# D is distance from particle (center of ribbon) to center of star
# L is full distance of ribbon projection

xvector = [x0]*n_particles
yvector = [y0]*n_particles
zvector = [z0]*n_particles

uvector = [0]*n_particles
vvector = [0]*n_particles
wvector = [0]*n_particles

xrotvec = [0]*n_particles
yrotvec = [0]*n_particles
zrotvec = [0]*n_particles

v_particle = 0.008 #actually is the v/V ratio, where v is v_particle and V is v_body
beta0 = 0.

betavec = [beta0]*n_particles

anglesrand = np.random.rand(n_particles, 2)

start = datetime.now()

for particle in xrange(n_particles):

    theta = math.acos(1. - 2.*anglesrand[particle, 0])
    phi = 2.*math.pi*anglesrand[particle, 1]

    u0 = v_particle*math.sin(theta)*math.cos(phi)+1./math.sqrt(2)
    v0 = v_particle*math.sin(theta)*math.sin(phi)+1./math.sqrt(2)
    w0 = v_particle*math.cos(theta)

    uvector[particle] = u0
    vvector[particle] = v0
    wvector[particle] = w0

dt = 0.1

#real units of time are acquired by using T * P / 2pi, where T is t in loop and P is period of 1 radian in orbit

timespan = np.linspace(0, 199.9, 2000) #end time = (dt*n_steps)-dt, or (0.1*1000)-0.1

for t in timespan:

    WDgrid = np.array([[0.]*201 for _ in range(201)])
    WDgrid.shape = (201, 201)

    for particle in xrange(n_particles):

        beta = betavec[particle]

        x = xvector[particle]
        y = yvector[particle]
        z = zvector[particle]

        u = uvector[particle]
        v = vvector[particle]
        w = wvector[particle]

        ix_0 = dt*(u)                                                         #dt*f(t, x, u)
        jx_0 = dt*(-x*(1-beta)/(math.sqrt(x**2+y**2+z**2))**3)                            #dt*g(t, x, u)
    
        iy_0 = dt*(v)                                                         #dt*f(t, x, u)
        jy_0 = dt*(-y*(1-beta)/(math.sqrt(x**2+y**2+z**2))**3)                                                #dt*g(t, x, u)

        iz_0 = dt*(w)                           
        jz_0 = dt*(-z*(1-beta)/(math.sqrt(x**2+y**2+z**2))**3)
        
        ix_1 = dt*(u + jx_0/2.)                                               #dt*f(t + dt/2., x + ix_0/2., u + jx_0/2.)
        jx_1 = dt*(-(x + ix_0/2.)*(1.-beta)/(math.sqrt((x + ix_0/2.)**2+(y + iy_0/2.)**2+(z + iz_0/2.)**2))**3)                        #dt*g(t + dt/2., x + ix_0/2., u + jx_0/2.)
    
        iy_1 = dt*(v + jy_0/2.)                                               #dt*f(t + dt/2., x + ix_0/2., u + jx_0/2.)
        jy_1 = dt*(-(y + iy_0/2.)*(1.-beta)/(math.sqrt((x + ix_0/2.)**2+(y + iy_0/2.)**2+(z + iz_0/2.)**2))**3)                        #dt*g(t + dt/2., x + ix_0/2., u + jx_0/2.)

        iz_1 = dt*(w + jz_0/2.)                                               
        jz_1 = dt*(-(z + iz_0/2.)*(1.-beta)/(math.sqrt((x + ix_0/2.)**2+(y + iy_0/2.)**2+(z + iz_0/2.)**2))**3)
        
        ix_2 = dt*(u + jx_1/2.)                                               #dt*f(t + dt/2., x + ix_1/2., u + jx_1/2.)
        jx_2 = dt*(-(x + ix_1/2.)*(1.-beta)/(math.sqrt((x + ix_1/2.)**2+(y + iy_1/2.)**2+(z + iz_1/2.)**2))**3)                        #dt*g(t + dt/2., x + ix_1/2., u + jx_1/2.)
    
        iy_2 = dt*(v + jy_1/2.)                                               #dt*f(t + dt/2., x + ix_1/2., u + jx_1/2.)
        jy_2 = dt*(-(y + iy_1/2.)*(1.-beta)/(math.sqrt((x + ix_1/2.)**2+(y + iy_1/2.)**2+(z + iz_1/2.)**2))**3)                        #dt*g(t + dt/2., x + ix_1/2., u + jx_1/2.)

        iz_2 = dt*(w + jz_1/2.)                                               #dt*f(
        jz_2 = dt*(-(z + iz_1/2.)*(1.-beta)/(math.sqrt((x + ix_1/2.)**2+(y + iy_1/2.)**2+(z + iz_1/2.)**2))**3)
        
        ix_3 = dt*(u + jx_2)                                                  #dt*f(t + dt, x + ix_2, u + jx_2)
        jx_3 = dt*(-(x + ix_2)*(1.-beta)/(math.sqrt((x + ix_2)**2+(y + iy_2)**2+(z + iz_2)**2))**3)                              #dt*g(t + dt, x + ix_2, u + jx_2)
    
        iy_3 = dt*(v + jy_2)                                                  #dt*f(t + dt, x + ix_2, u + jx_2)
        jy_3 = dt*(-(y + iy_2)*(1.-beta)/(math.sqrt((x + ix_2)**2+(y + iy_2)**2+(z + iz_2)**2))**3)                              #dt*g(t + dt, x + ix_2, u + jx_2)

        iz_3 = dt*(w + jz_2)                                               
        jz_3 = dt*(-(z + iz_2)*(1.-beta)/(math.sqrt((x + ix_2)**2+(y + iy_2)**2+(z + iz_2)**2))**3)
    
        x = x + (1./6.)*(ix_0 + 2*ix_1 + 2*ix_2 + ix_3)
        u = u + (1./6.)*(jx_0 + 2*jx_1 + 2*jx_2 + jx_3)
    
        y = y + (1./6.)*(iy_0 + 2*iy_1 + 2*iy_2 + iy_3)
        v = v + (1./6.)*(jy_0 + 2*jy_1 + 2*jy_2 + jy_3)

        z = z + (1./6.)*(iz_0 + 2*iz_1 + 2*iz_2 + iz_3)
        w = w + (1./6.)*(jz_0 + 2*jz_1 + 2*jz_2 + jz_3)

        xvector[particle] = x
        yvector[particle] = y
        zvector[particle] = z

        uvector[particle] = u
        vvector[particle] = v
        wvector[particle] = w

print datetime.now()-start
