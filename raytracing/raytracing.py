#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as LA

'''
Code from Jevon Longdell for ray tracing
based on G. H. Spencer and M. V. R. K. Murty "General Ray-Tracing Procedure" Journal of the Optical Society ofAMeria 52(6) 1962.
'''

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