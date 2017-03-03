import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np

#x'' = -GMx/sqrt(x^2+y^2)^3
#y'' = -GMy/sqrt(x^2+y^2)^3
#
#x' = u                     1st eq, f(t, x, y)
#u' = -GMx/(x^2+y^2)^3      2nd eq, g(t, x, y)
#
#y' = v
#v' = -GMy/(x^2+y^2)^3


n_particles = 1000  #4.5 minutes for 100 particles, also 4.5 minutes for 1000 particles

t = 0.0

x0 = 1./math.sqrt(2)
y0 = -1./math.sqrt(2)
z0 = 0.0

# impact parameter, b, is the z-coordinate
# D is distance from particle (center of ribbon) to center of star
# L is full distance of ribbon projection

lightcurve = []

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
beta0 = 0.0001

betavec = [beta0]*n_particles

anglesrand = np.random.rand(n_particles, 2)
betaran = np.random.rand(n_particles)

for particle in xrange(n_particles):

    #temptheta = math.acos(1 - 0.0038*anglesrand[particle, 0]) #0.0038 is (1 - cosine of 5 degrees), the radius of the cone at the top and bottom of the z-axis

    #if anglesrand[particle, 0] < 0.5:
    #    theta = temptheta
    #if anglesrand[particle, 0] >= 0.5:
    #    theta = math.pi - temptheta

    theta = math.acos(1. - 2.*anglesrand[particle, 0])
    #theta = 0           #testing if at theta = 0 there should be little to no expansion in x or y, only z-direction
    phi = 2.*math.pi*anglesrand[particle, 1]

    u0 = v_particle*math.sin(theta)*math.cos(phi)+1./math.sqrt(2)
    v0 = v_particle*math.sin(theta)*math.sin(phi)+1./math.sqrt(2)
    w0 = v_particle*math.cos(theta)

    uvector[particle] = u0
    vvector[particle] = v0
    wvector[particle] = w0

    #betavec[particle] += 0.19*betaran[particle]

    #if betaran[particle] < 0.5:
    #    betavec[particle] = 0.1

dt = 0.1

v_ratiostr = "{} = {}".format("v/V", v_particle)
#betastr = "{} = {}".format("Beta", beta)

#real units of time are acquired by using T * P / 2pi, where T is t in loop and P is period of 1 radian in orbit

timespan = np.linspace(0, 99.9, 1000) #end time = (dt*n_steps)-dt, or (0.1*1000)-0.1

for t in timespan:

    count = 0

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

        sine = math.sin(t)
        cosine = math.cos(t)

        xvector[particle] = x
        yvector[particle] = y
        zvector[particle] = z

        uvector[particle] = u
        vvector[particle] = v
        wvector[particle] = w

        xrotvec[particle] = x*cosine + y*sine
        yrotvec[particle] = y*cosine - x*sine
        zrotvec[particle] = z

        if ((x**2+z**2) < (0.011)**2 and y < 0.): 
            count += 1 

    lightcurve.append(count)

    circle1 = plt.Circle((0,0), 0.011, color='cyan')
    circle2 = plt.Circle((0,0), 0.011, color='cyan')
    circle3 = plt.Circle((0,0), 0.011, color='cyan')
    circle4 = plt.Circle((0,0), 0.011, color='cyan')
    circle5 = plt.Circle((0,0), 1.0, color='red', fill=False)
    circle6 = plt.Circle((0,0), 1.0, color='red', fill=False)
    circle7 = plt.Circle((0,0), 1.0, color='red', fill=False)
    circle8 = plt.Circle((0,0), 1.0, color='red', fill=False)

    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim([-1.2*1.2, 1.2*1.2])
    ax.set_ylim([-1.2, 1.2])
    ax.add_artist(circle1)
    ax.add_artist(circle5)
    plt.plot(xrotvec, yrotvec, 'ro', ms=1)
    plt.figtext(0.705, 0.825, v_ratiostr, fontsize = 'large')
    #plt.figtext(0.705, 0.785, betastr, fontsize = 'large')
    figname = '{}{}'.format('RKExplodingOrbit/xyimgrot', int(t*10))
    plt.savefig(figname)
    plt.close()

    plt.axis('equal')
    ax3 = plt.gca()
    ax3.set_xlim([-1.2*1.2, 1.2*1.2])
    ax3.set_ylim([-1.2, 1.2])
    ax3.add_artist(circle2)
    ax3.add_artist(circle6)
    plt.plot(xvector, yvector, 'ro', ms=1)
    plt.figtext(0.705, 0.825, v_ratiostr, fontsize = 'large')
    #plt.figtext(0.705, 0.785, betastr, fontsize = 'large')
    figname3 = '{}{}'.format('RKExplodingOrbit/xyimg', int(t*10))
    plt.savefig(figname3)
    plt.close()

    plt.axis('equal')
    ax2 = plt.gca()
    ax2.set_xlim([-1.2*1.2, 1.2*1.2])
    ax2.set_ylim([-1.2, 1.2])
    ax2.add_artist(circle3)
    ax2.add_artist(circle7)
    plt.plot(xrotvec, zrotvec, 'ro', ms=1)
    plt.figtext(0.705, 0.825, v_ratiostr, fontsize = 'large')
    #plt.figtext(0.705, 0.785, betastr, fontsize = 'large')
    figname2 = '{}{}'.format('RKExplodingOrbit/xzimgrot', int(t*10))
    plt.savefig(figname2)
    plt.close()

    plt.axis('equal')
    ax4 = plt.gca()
    ax4.set_xlim([-1.2*1.2, 1.2*1.2])
    ax4.set_ylim([-1.2, 1.2])
    ax4.add_artist(circle4)
    ax4.add_artist(circle8)
    plt.plot(xvector, zvector, 'ro', ms=1)
    plt.figtext(0.705, 0.825, v_ratiostr, fontsize = 'large')
    #plt.figtext(0.705, 0.785, betastr, fontsize = 'large')
    figname4 = '{}{}'.format('RKExplodingOrbit/xzimg', int(t*10))
    plt.savefig(figname4)
    plt.close()

    t = t + dt

#period = 2*math.pi

#foldarray = [0] * 100
#exposure = [0] * 100
#for i in range(0, len(lightcurve)):
#    frac = (lightcurve[i]%period)/period
#    foldarray[int(frac*100)] += lightcurve[i]
#    exposure[int(frac*100)] += 1
#foldedcurve = np.divide(foldarray, exposure)
foldedcurvebins = np.linspace(0, 100, len(lightcurve))
norm = [1000]*len(lightcurve)
for i in xrange(len(lightcurve)):
    lightcurve[i] = norm[i] - lightcurve[i]

plt.step(foldedcurvebins, lightcurve, 'k')
plt.savefig('LightCurve')
plt.close()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#    iy_0 = dt*(v)                                                         #dt*f(t, x, u)
#    jy_0 = dt*(-k*y - b*v)                                                #dt*g(t, x, u)
#    iy_1 = dt*(v + jy_0/2.)                                               #dt*f(t + dt/2., x + ix_0/2., u + jx_0/2.)
#    jy_1 = dt*(-k*(y + iy_0/2.) - b*(v + jy_0/2.))                        #dt*g(t + dt/2., x + ix_0/2., u + jx_0/2.)
#    iy_2 = dt*(v + jy_1/2.)                                               #dt*f(t + dt/2., x + ix_1/2., u + jx_1/2.)
#    jy_2 = dt*(-k*(y + iy_1/2.) - b*(v + jy_1/2.))                        #dt*g(t + dt/2., x + ix_1/2., u + jx_1/2.)
#    iy_3 = dt*(v + jy_2)                                                  #dt*f(t + dt, x + ix_2, u + jx_2)
#    jy_3 = dt*(-k*(y + iy_2) - b*(v + jy_2))                              #dt*g(t + dt, x + ix_2, u + jx_2)
#    y = y + (1/6.)*(iy_0 + 2*iy_1 + 2*iy_2 + iy_3)
#    v = v + (1/6.)*(jy_0 + 2*jy_1 + 2*jy_2 + jy_3)

#x'' = -(2*a*x + b*y) - d(x' + y')
#y'' = -(b*x + 2*c*y) - d(x' + y')

#x' = u                              1st eq, f(t, x)
#u' = -(2*a*x + b*y) - d(u + v)      2nd eq, g(t, x)

#y' = v
#v' = -(b*x + 2*c*y) - d(u + v)

#m = 5.0
#t = 1.0
#x = 0.0
#y = 0.0
#u = 0.5
#v = 0.5
##k = 0.5/m
#drag = 0.15/m #drag
#a = 1.0/m #general constants
#b = 0.5/m
#c = 2.0/m
#
#dt = 0.05
#
#tvector = []
#xvector = []
#yvector = []
#
#tvector.append(t)
#xvector.append(x)
#yvector.append(y)
#
#while (t < 200.0):
#
#    ix_0 = dt*(u)                                 
#    jx_0 = dt*(-(2*a*x + b*y) - drag*(u + v))    
#
#    iy_0 = dt*(v)                                                         
#    jy_0 = dt*(-(b*x + 2*c*y) - drag*(u + v))                                                                   
#    
#    ix_1 = dt*(u + jx_0/2.)                       
#    jx_1 = dt*(-(2*a*(x + ix_0/2.) + b*(y + iy_0/2.)) - drag*((u + jx_0/2.) + (v + jy_0/2.)))
#    
#    iy_1 = dt*(v + jy_0/2.)                                               
#    jy_1 = dt*(-(b*(x + ix_0/2.) + 2*c*(y + iy_0/2.)) - drag*((u + jx_0/2.) + (v + jy_0/2.)))                        
#
#    ix_2 = dt*(u + jx_1/2.)                       
#    jx_2 = dt*(-(2*a*(x + ix_1/2.) + b*(y + iy_1/2.)) - drag*((u + jx_1/2.) + (v + jy_1/2.)))
#
#    iy_2 = dt*(v + jy_1/2.)                                               
#    jy_2 = dt*(-(b*(x + ix_1/2.) + 2*c*(y + iy_1/2.)) - drag*((u + jx_1/2.) + (v + jy_1/2.)))                       
#    
#    ix_3 = dt*(u + jx_2)                          
#    jx_3 = dt*(-(2*a*(x + ix_2) + b*(y + iy_2)) - drag*((u + jx_2) + (v + jy_2))) 
#
#    iy_3 = dt*(v + jy_2)                                                  
#    jy_3 = dt*(-(b*(x + ix_2) + 2*c*(y + iy_2)) - drag*((u + jx_2) + (v + jy_2)))                            
#
#    t = t + dt
#    x = x + (1/6.)*(ix_0 + 2*ix_1 + 2*ix_2 + ix_3)
#    u = u + (1/6.)*(jx_0 + 2*jx_1 + 2*jx_2 + jx_3)
#
#    y = y + (1/6.)*(iy_0 + 2*iy_1 + 2*iy_2 + iy_3)
#    v = v + (1/6.)*(jy_0 + 2*jy_1 + 2*jy_2 + jy_3)
#    
#    
#    tvector.append(t)
#    xvector.append(x)
#    yvector.append(y)
#
##plt.plot(tvector, xvector, 'k')
##plt.show()
#
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.plot(xvector, yvector, tvector, label='parametric curve')
#plt.xlabel('X-axis')
#plt.ylabel('Y-axis')
#plt.show()