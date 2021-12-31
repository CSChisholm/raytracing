#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as LA

#Code from Jevon Longdell for ray tracing

class Ray:
    '''A ray object to pass to trace()'''
    def __init__(self):
        self.direction = np.array([0,0,1]) #The propagation direction of the ray
        self.r0 = np.array([0,0,0]) #The start point of the ray
        self.length = 0 #The length of the ray
        self.wavelength = 589 #Wavelength in nm
        
    def plot(self):
        '''Plot the ray'''
        r1 = self.r0 + self.length*self.direction
        plt.plot((self.r0[2],r1[2]),(self.r0[0],r1[0]),'b')

class SphericalSurface:
    '''Spherical surface'''
    def __init__(self):
        self.Z0 = 0 #Surface position along optical axis
        self.curv = 0 #Curvature of surface
        self.n1 = 1 #Refractive index to lef tof surface
        self.n2 = 2 #Refractive index to right of surface
        self.app = 6 #Apperature of surface
        self.ccurv = np.array([0,0,self.Z0+1/self.curv]) #Centre of curvature of the surface

    def Ffunc(self,r):
        '''This is the F of eq. 1 in Spencer and Murty'''
        F = np.dot(r-self.ccurv,r-self.ccurv)-(1/self.curv)**2 #
        return F
    
    def gradFfunc(self,r):
        return 2*(r - self.ccurv)
    
    def Zfunc(self,ssq):
        #not used in calcs
        c = self.curv
        Z = self.Z0 + c*ssq/(1+np.sqrt(1-c*c*ssq))
        return Z

class Surface:
    
    def __init__(self):
        self.Z0 = 0 #Surface position along optical axis
        self.curv = 0 #Curvature of surface
        self.kappa = 1 #=1+k, k = Conic coefficient
        self.Aparams = np.array([]) #Polynomial coefficients of asphere
        self.n1 = 1.0 #Refractive index to lef tof surface
        self.n2 = 1.0 #Refractive index to right of surface
        self.app = 6 #Apperature of surface
        
    def Zfunc(self,ssq):
        c = self.curv
        Z = self.Z0 + c*ssq/(1+np.sqrt(1-self.kappa*c*c*ssq))
        for k in range(len(self.Aparams)):
            Z = Z+self.Aparams[k]*ssq**(2+k)
        return Z

    def Ffunc(self,r):
        '''This is the F of eq. 1 in Spencer and Murty'''
        x = r[0]
        y = r[1]
        z = r[2]
        ssq = x*x+y*y
        F = z-self.Zfunc(ssq)
        return F

    def gradFfunc(self,r):
        x = r[0]
        y = r[1]
        ssq = x*x+y*y
        c = self.curv
        E = c/np.sqrt(1-self.kappa*c*c*ssq)
        for k in range(len(self.Aparams)):
            E = E + 2*(2+k)*self.Aparams[k]*ssq**(1+k) #this could be wrong
        Fx = -x*E
        Fy = -y*E
        Fz = 1
        return np.array([Fx,Fy,Fz])



def trace(ray_bundle,surf):
    '''Perform the ray tracing procedure'''
    new_ray_bundle = []

    #calculates the center of curvature of the surface
    if abs(surf.curv)>1/100:
        ccurv = np.array([0,0,surf.Z0+1./surf.curv])
        print(f'{ccurv}')
    else:
        print('nope\n\n')
        ccurv = None

    count=0
    for ray in ray_bundle:
    #first determine value for s where the ray crosses the z = z0 plane
        #ccurv for the lens
        count+=1
        print(f'Ray = {count}')
        if (ccurv is not None):
            tv = ray.r0 - ccurv #temp vector
            b = 2*np.dot(ray.direction,tv)
            c = np.dot(tv,tv) - 1/surf.curv**2
            s1 = (-b - np.sqrt(b**2 - 4*c))/2
            s2 = (-b + np.sqrt(b**2 - 4*c))/2
            if (surf.curv>0):
                s = s1
            else:
                s = s2
            print([s1, s2, s,b**2-4*c])
            if ((b**2 - 4*c)<0):
                s = (ray.r0[2]-surf.Z0)/ray.direction[2]
        else:
            s = (ray.r0[2] - surf.Z0)/ray.direction[2]
        for k in range(30):
            r = ray.r0 + s*ray.direction
            dFds = np.dot(ray.direction,surf.gradFfunc(r))
            snew = s - surf.Ffunc(r)/dFds
            if (np.absolute(s-snew)<1e-9):
                break
            s = snew
        if k>25:
            #here we assume that the ray has missed the surface so for the new ray make r0=r0 and length of old ray zero
            ray.length = 0
            newray = Ray()
            newray.r0 = ray.r0
            newray.direction = ray.direction
        else:
            ray.length = snew
            r = ray.r0 + snew*ray.direction
            normal = surf.gradFfunc(r)
            normal = normal/LA.norm(normal)

            newray = Ray()
            newray.r0 = r
            print(f'{np.dot(ray.direction,normal)}')
            a = surf.n1/surf.n2*np.dot(ray.direction,normal)
            b = ((surf.n1/surf.n2)**2-1)
            gamma = -b/(2*a)
            for k in range(30):
                gammanew = (gamma**2-b)/(2*(gamma+a))
                if (np.absolute(gamma-gammanew)<1e-8):
                    gamma = gammanew
                    break
                gamma = gammanew
            assert(k<25)
            newray.direction = surf.n1/surf.n2*ray.direction+gamma*normal
        assert(np.absolute(np.absolute(LA.norm(newray.direction)-1))<1e-8)
        new_ray_bundle.append(newray)
    return new_ray_bundle

def partialsum(vector,ind):
    '''A function for partial sum'''
    return np.sum(vector[:(ind+1)])

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
        r = Ray()
        r.r0 = np.array([rayoffset,0,0])
        r.direction = np.array([np.sin(x),0,np.cos(x)])
        rays.append(r)
    
    #Set surfaces
    offset = 0
    surfs = []
    
    offset+=atomposition
    s1 = Surface()
    s1.Z0 = offset
    s1.n2 = nBK7
    s1.app = asphapp
    surfs.append(s1)
    
    offset+=asphcenthick
    s2 = Surface()
    s2.curv = -1/asphR
    s2.kappa = 1+asphk
    s2.n1 = nBK7
    s2.app = asphapp
    s2.Aparams = asphparams
    s2.Z0 = offset
    surfs.append(s2)
    
    offset+=ww
    s3 = Surface()
    s3.curv = 1/surfrad2
    s3.Z0 = offset
    s3.app = app2
    s3.n2 = nBK7
    surfs.append(s3)
    
    offset+=centhick2
    s4 = Surface()
    s4.Z0 = offset
    s4.n1 = nBK7
    s4.app = app2
    surfs.append(s4)
    
    if (inclen3==1):
        offset+=ww2
        if (len3ori==1):
            s5 = Surface()
            s5.curv = 1/surfrad3
            s5.Z0 = offset
            s5.n2 = nBK7
            s5.app = app3
            surfs.append(s5)
            
            offset+=centhick3
            s6 = Surface()
            s6.Z0 = offset
            s6.n1 = nBK7
            s6.app = app3
            surfs.append(s6)
        else:
            s5 = Surface()
            s5.Z0 = offset
            s5.n2 = nBK7
            s5.app = app3
            surfs.append(s5)
            
            offset+=centhick3
            s6 = Surface()
            s6.curv = -1/surfrad3
            s6.Z0 = offset
            s6.n1 = nBK7
            s6.app = app3
            surfs.append(s6)
    
    send = Surface()
    send.Z0 = camerapos
    send.app = 512*24e-3
    surfs.append(send)
    
    #Calculations
    plt.figure()
    for s in surfs:
        xvals = np.linspace(-s.app,s.app,100)/2
        zvals = 0*xvals;
        for k in range(len(xvals)):
            zvals[k] = s.Zfunc(xvals[k]**2)
        plt.plot(zvals,xvals,'r')
    
        newrays = trace(rays,s)
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