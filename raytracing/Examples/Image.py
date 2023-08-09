#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from raytracing import raytracing as rt

def lensparam(fin):
    '''Parameters of 2" plano-convex lenses from Thorlabs'''
    if fin==75:
        surfradout = 38.6 #Radius of curvature in mm
        centhickout = 12.5 #Centre thickness in mm
    elif fin==100:
        surfradout = 51.5 #Radius of curvature in mm
        centhickout = 9.7 #Centre thickness in mm
    elif fin==200:
        surfradout = 103.0 #Radius of curvature in mm
        centhickout = 6.2 #Centre thickness in mm
    elif fin==300:
        surfradout = 154.5 #Radius of curvature in mm
        centhickout = 5.1 #Centre thickness in mm
    elif fin==400:
        surfradout = 206.0 #Radius of curvature in mm
        centhickout = 4.6 #Centre thickness in mm
    elif fin==500:
        surfradout = 257.5 #Radius of curvature in mm
        centhickout = 4.3 #Centre thickness in mm
    elif fin==1000:
        surfradout = 515.1 #Radius of curvature in mm
        centhickout = 3.6 #Centre thickness in mm
    else:
        surfradout = float('inf') #Radius of curvature in mm
        centhickout = 5.1 #Centre thickness in mm
    return surfradout, centhickout

def main():
    plt.close('all')
    
    #Some parameters
    nBK7 = 1.517
    
    #AL2550H parameters
    asphR = 25.56 #Surface curvature in mm
    asphk = -1.01 #Conic coefficient
    asphparams = np.array([3.270398e-6,7.7205335e-10,1.6304727e-13])
    asphcenthick = 6.0 #Centre thickness of asphere in mm
    asphapp = 25.0 #Aperture of asphere in mm
    
    #Lens2
    f2 = 300 #Approximate focal length in mm
    surfrad2, centhick2 = lensparam(f2)
    app2 = 50.8 #Aperture of second lens in mm
    
    #lens3
    inclen3 = 0 #1 = include, anything else = don't include
    len3ori = 0 #1 = curve towards atoms, else = plane toward atoms
    f3 = 75
    surfrad3, centhick3 = lensparam(f3)
    app3 = 50.8
    
    #Other parameters
    atomposition = 43 #Distance from atoms to first lens
    ww = 440 #Distance betwqeen lenses in mm
    ww2 = 220 #Distance from second to third lens in mm
    camerapos = 1100 #Distance from atoms to camera in mm
    
    #Set source rays
    maxang = 0.1
    numrays = 22
    rays = []
    
    rayoffset = 0.2*0
    for x in np.linspace(-maxang,maxang,numrays):
        r = rt.Ray()
        r.r0 = np.array([rayoffset,0,0])
        r.direction = np.array([np.sin(x),0,np.cos(x)])
        rays.append(r)
    
    #Set surfaces
    offset = 0
    surfs = []
    
    offset+=atomposition
    s1 = rt.Surface()
    s1.Z0 = offset
    s1.n2 = nBK7
    s1.app = asphapp
    surfs.append(s1)
    
    offset+=asphcenthick
    s2 = rt.Surface()
    s2.curv = -1/asphR
    s2.kappa = 1+asphk
    s2.n1 = nBK7
    s2.app = asphapp
    s2.Aparams = asphparams
    s2.Z0 = offset
    surfs.append(s2)
    
    offset+=ww
    s3 = rt.Surface()
    s3.curv = 1/surfrad2
    s3.Z0 = offset
    s3.app = app2
    s3.n2 = nBK7
    surfs.append(s3)
    
    offset+=centhick2
    s4 = rt.Surface()
    s4.Z0 = offset
    s4.n1 = nBK7
    s4.app = app2
    surfs.append(s4)
    
    if (inclen3==1):
        offset+=ww2
        if (len3ori==1):
            s5 = rt.Surface()
            s5.curv = 1/surfrad3
            s5.Z0 = offset
            s5.n2 = nBK7
            s5.app = app3
            surfs.append(s5)
            
            offset+=centhick3
            s6 = rt.Surface()
            s6.Z0 = offset
            s6.n1 = nBK7
            s6.app = app3
            surfs.append(s6)
        else:
            s5 = rt.Surface()
            s5.Z0 = offset
            s5.n2 = nBK7
            s5.app = app3
            surfs.append(s5)
            
            offset+=centhick3
            s6 = rt.Surface()
            s6.curv = -1/surfrad3
            s6.Z0 = offset
            s6.n1 = nBK7
            s6.app = app3
            surfs.append(s6)
    
    send = rt.Surface()
    send.Z0 = camerapos
    send.app = 512*24e-3
    surfs.append(send)
    
    #Calculations
    plt.figure()
    for s in surfs:
        xvals = np.linspace(-s.app,s.app,100)/2
        zvals = [s.Zfunc(xval**2) for xval in xvals]
        plt.plot(zvals,xvals,'C3')
        newrays = rt.trace(rays,s)
        for r in rays:
            r.plot()
        rays=newrays
    for r in rays:
        r.plot()
    plt.axis([0,camerapos+50,-50,50])
    plt.xlabel('z [mm]')
    plt.ylabel('x, y [mm]')
    plt.show()

if (__name__=='__main__'):
    main()