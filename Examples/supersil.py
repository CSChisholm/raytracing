#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from raytracing import raytracing as rt

def main():
    plt.close('all')
    rays = []
    
    for x in np.linspace(1.5,-1.5,20):
        r = rt.Ray()
        r.r0 = np.array([x,0,-8.5])
        rays.append(r)
    #set surfaces
    surfs = []
    
    #parameters of the lightpath 352330 aspheric
    as1curv = 1/-3.20561
    as1kappa = 1+-12.418013
    as1Aparams = np.array([9.00531e-3,-1.359752e-3,1.136638e-4,-4.278925e-6])
    athickness = 3.19
    awd = 1.76
    as2curv = 1/2.74797
    as2kappa = 1+-0.542698
    as2Aparams = np.array([-3.19546e-4,-4.397785e-5,1.842256e-5,-1.566446e-6])
    
    nglass = 1.595
    nyso = 1.8
    ballradius = 0.625
    
    zshort = ballradius/nyso
    zlong = ballradius*nyso
    
    offset = +zlong-zshort
    
    s1 = rt.Surface()
    s1.curv = as2curv
    s1.kappa = as2kappa
    s1.Aparams = as2Aparams
    s1.Z0 = -1.76-3.19+offset
    s1.n2 = nglass
    
    s2 = rt.Surface()
    s2.curv = as1curv
    s2.kappa = as1kappa
    s2.Aparams = as1Aparams
    s2.Z0 = -1.76+offset
    s2.n1 = nglass
    
    s3 = rt.SphericalSurface()
    s3.curv = 1/ballradius
    s3.Z0 = -ballradius-zshort
    s3.n2 = nyso
    s3.app = 2*ballradius
    
    s4 = rt.Surface()
    s4.Z0 = 0.01
    
    surfs.append(s1)
    surfs.append(s2)
    surfs.append(s3)
    surfs.append(s4)
    
    plt.figure()
    for s in surfs:
        xvals = np.linspace(-s.app,s.app,100)/2
        zvals = [s.Zfunc(xval**2) for xval in xvals]
        plt.plot(zvals,xvals,'r')
    
        newrays = rt.trace(rays,s)
        for r in rays:
            r.plot()
           
        rays=newrays
    
    for r in rays:
        r.plot()
    
    plt.axis([-5.5,5.5,-3,3])    
    plt.axis('equal')
    plt.savefig('supersil.pdf')
    plt.show()

if (__name__=='__main__'):
    main()