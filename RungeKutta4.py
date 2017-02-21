import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

#x'' = -GMx/sqrt(x^2+y^2)^3
#y'' = -GMy/sqrt(x^2+y^2)^3
#
#x' = u                     1st eq, f(t, x, y)
#u' = -GMx/(x^2+y^2)^3      2nd eq, g(t, x, y)
#
#y' = v
#v' = -GMy/(x^2+y^2)^3


#m = 5.0
t = 1.0
x = 1.0
y = 0.0
u = 0.0
v = 1.3
#k = 0.5/m
#b = 0.2/m
G = math.pow(6.67, -8)
M = math.pow(1.99, 33)

omega = math.sqrt(G*M/(x**2+y**2))

dt = 0.1

#real units of time are acquired by using T * P / 2pi, where T is t in loop and P is period of 1 radian in orbit

tvector = []
xvector = []
yvector = []
tvector.append(t)
xvector.append(x)
yvector.append(y)

while (t < 200.0):

    ix_0 = dt*(u)                                                         #dt*f(t, x, u)
    jx_0 = dt*(-x/(math.sqrt(x**2+y**2))**3)                            #dt*g(t, x, u)

    iy_0 = dt*(v)                                                         #dt*f(t, x, u)
    jy_0 = dt*(-y/(math.sqrt(x**2+y**2))**3)                                                #dt*g(t, x, u)
    
    ix_1 = dt*(u + jx_0/2.)                                               #dt*f(t + dt/2., x + ix_0/2., u + jx_0/2.)
    jx_1 = dt*(-(x + ix_0/2.)/(math.sqrt((x + ix_0/2.)**2+(y + iy_0/2.)**2))**3)                        #dt*g(t + dt/2., x + ix_0/2., u + jx_0/2.)

    iy_1 = dt*(v + jy_0/2.)                                               #dt*f(t + dt/2., x + ix_0/2., u + jx_0/2.)
    jy_1 = dt*(-(y + iy_0/2.)/(math.sqrt((x + ix_0/2.)**2+(y + iy_0/2.)**2))**3)                        #dt*g(t + dt/2., x + ix_0/2., u + jx_0/2.)
    
    ix_2 = dt*(u + jx_1/2.)                                               #dt*f(t + dt/2., x + ix_1/2., u + jx_1/2.)
    jx_2 = dt*(-(x + ix_1/2.)/(math.sqrt((x + ix_1/2.)**2+(y + iy_1/2.)**2))**3)                        #dt*g(t + dt/2., x + ix_1/2., u + jx_1/2.)

    iy_2 = dt*(v + jy_1/2.)                                               #dt*f(t + dt/2., x + ix_1/2., u + jx_1/2.)
    jy_2 = dt*(-(y + iy_1/2.)/(math.sqrt((x + ix_1/2.)**2+(y + iy_1/2.)**2))**3)                        #dt*g(t + dt/2., x + ix_1/2., u + jx_1/2.)
    
    ix_3 = dt*(u + jx_2)                                                  #dt*f(t + dt, x + ix_2, u + jx_2)
    jx_3 = dt*(-(x + ix_2)/(math.sqrt((x + ix_2)**2+(y + iy_2)**2))**3)                              #dt*g(t + dt, x + ix_2, u + jx_2)

    iy_3 = dt*(v + jy_2)                                                  #dt*f(t + dt, x + ix_2, u + jx_2)
    jy_3 = dt*(-(y + iy_2)/(math.sqrt((x + ix_2)**2+(y + iy_2)**2))**3)                              #dt*g(t + dt, x + ix_2, u + jx_2)

    t = t + dt

    x = x + (1/6.)*(ix_0 + 2*ix_1 + 2*ix_2 + ix_3)
    u = u + (1/6.)*(jx_0 + 2*jx_1 + 2*jx_2 + jx_3)

    y = y + (1/6.)*(iy_0 + 2*iy_1 + 2*iy_2 + iy_3)
    v = v + (1/6.)*(jy_0 + 2*jy_1 + 2*jy_2 + jy_3)
    
    
    tvector.append(t)
    xvector.append(x)
    yvector.append(y)

plt.axis('equal')
plt.plot(xvector, yvector, 'k')
plt.show()
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.plot(xvector, yvector, tvector, label='parametric curve')
#plt.show()

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