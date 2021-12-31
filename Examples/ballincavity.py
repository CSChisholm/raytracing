#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from raytracing import raytracing as rt

def main():
    #################################################
    ## setup the problem
    #######################################
    
    plt.close('all')
    plt.ion()
    
    b = rt.Ball()
    b.n_inside=1.6
    b.n_outside=1
    b.label = 'centralball'
    
    m1 = rt.Ball()
    m1.radius = 20
    m1.n_inside=1
    m1.n_outside = np.inf #it's a mirror
    m1.center = np.array([0,0,10])
    m1.label = 'mirror1'
    
    m2 = rt.Ball()
    m2.radius = 20
    m2.center = np.array([0,0,-10])
    m2.n_inside=1
    m2.n_outside = np.inf #it's a mirror
    m2.label = 'mirror2'
    plt.figure()
    b.plot()
    m1.plot()
    m2.plot()
    
    rays = []
    for th in np.linspace(-0.1,0.1,5):
        r = rt.Ray()
        r.r0 = np.array([0,0,-3])
        r.direction = np.array([np.sin(th),0,np.cos(th)])
        r.length = 1
        rays.append(r)
    
    surfs = [b,m1,m2]
    plt.figure()
    for k in range(4):
        newrays = rt.traceball(rays,surfs)
        for r in rays:
            r.plot()
        rays=newrays    
    plt.axis((-11,11,-3,3))
    plt.show()
     
if (__name__=='__main__'):
    main()