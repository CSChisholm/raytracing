#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as LA
from scipy.optimize import brentq

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
        self.curv = 1 #Curvature of surface
        self.n1 = 1 #Refractive index to lef tof surface
        self.n2 = 2 #Refractive index to right of surface
        self.app = 6 #Apperature of surface
    
    def ccurv(self):
        '''Centre of curvature of the surface'''
        return np.array([0,0,self.Z0+1/self.curv])

    def Ffunc(self,r):
        '''This is the F of eq. 1 in Spencer and Murty'''
        F = np.dot(r-self.ccurv(),r-self.ccurv())-(1/self.curv)**2 #
        return F
    
    def gradFfunc(self,r):
        return 2*(r - self.ccurv())
    
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
    
    def ccurv(self):
        '''Centre of curvature of the surface'''
        return np.array([0,0,self.Z0+1/self.curv])    
    
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

class Ball:
    #the inside is defined as the region F>0
    def __init__(self):
        self.center = np.array([0,0,0])
        self.radius = 1
        self.n_inside = 1.0
        self.n_outside = 1.8
        self.label = 'round thing'
    
    def __str__(self):
        return f'Ball with label {self.label}'
    
    def plot(self):
        thetavals = np.linspace(0,2*np.pi,1024)
        zvals = self.center[2]+self.radius*np.sin(thetavals)
        xvals = self.center[0]+self.radius*np.cos(thetavals)
        plt.plot(zvals,xvals,'r')
    
    def Ffunc(self,r):
        '''returns radius^2 - |r-center|^2'''
        return self.radius**2 - np.dot(r-self.center,r-self.center)
    
    def gradFfunc(self,r):
        return 2*(self.center-r)
                
    def find_intersect(self,ray):
        '''does the ray intercept the surface if so how far along ray? if not return -1'''
        s0 = 1e-5 #start this distance along the ray so that we don't get confused by the surface we have just passed thru
        if (self.Ffunc(ray.r0+s0*ray.direction)>0): #do we start of inside the sphere
            #code for if we do:
            #if you go three radii in any direction you are definatley outside the sphere so this point is definately outside the sphere
            soutside = 3*self.radius
            #now we have values of s that bracket the intercept we can find it
            s = brentq(lambda x: self.Ffunc(ray.r0+x*ray.direction),s0,soutside)
            return s
        else:
            #code for if we start outside sphere
            #s for closest point to the center
            s_closest = np.dot(self.center-ray.r0,ray.direction)
            if s_closest < s0:
                return -1 #we are outside and heading away from center
            #f value for closest point to the center
            F_closest = self.Ffunc(ray.r0+s_closest*ray.direction)
            if (F_closest<0): #ray never hits sphere
                return -1
            else:
                #s_closest and 0 are on different sides of the sphere
                s = brentq(lambda x: self.Ffunc(ray.r0+x*ray.direction),0,s_closest)
                return s
        #we should never reach this point
        #this is a way to make sure we know if we do
        assert(False) 

def trace(ray_bundle,surf):
    '''Perform the ray tracing procedure'''
    new_ray_bundle = []

    #calculates the center of curvature of the surface
    if (abs(surf.curv)>1/100):
        ccurv = surf.ccurv()
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

def traceball(ray_bundle,surfs):
    new_ray_bundle=[]

    for ray in ray_bundle:
        # work out the next surface for the ray to hit
        # do this by working out all the surfaces the ray
        # intesects and choosing the smallest value of s
        smin=np.inf
        nextsurf = None
        for surf in surfs:
            s = surf.find_intersect(ray)
            #print(surf,s)
            if s>=0:
                if s<smin:
                    smin=s
                    nextsurf=surf
        if not(np.isfinite(smin)):
            print("Ray missed all the surfaces")
            assert(False)
        ray.length=smin
        r = ray.r0+(smin)*ray.direction
        normal = nextsurf.gradFfunc(r)
        normal = normal/LA.norm(normal)
        newray = Ray()
        newray.r0 = r
        newray.length=0.2
        #this chooses the direction of the new ray
        n1 = nextsurf.n_outside
        n2 = nextsurf.n_inside
        # if either of the refractive indicies are infinite
        # then act like is't s a mirror
        # otherwise use snells law
        if (np.isfinite(nextsurf.n_outside) and np.isfinite(nextsurf.n_inside)):
            # donesn't look like snells law
            # but see the a paper
            if (np.dot(ray.direction,normal)<0):
                n2 = nextsurf.n_outside
                n1 = nextsurf.n_inside
            else:
                n1 = nextsurf.n_outside
                n2 = nextsurf.n_inside     
            a = n1/n2*np.dot(ray.direction,normal)
            b = ((n1/n2)**2-1)
            gamma = -b/(2*a)
            for k in range(30):
                gammanew = (gamma**2-b)/(2*(gamma+a))
                if abs(gamma-gammanew) <1e-8:
                    gamma = gammanew
                    break
                gamma=gammanew
            #make sure iterative approach has converged 
            assert(k<25)
            newray.direction = n1/n2*ray.direction+gamma*normal
            assert(np.absolute(LA.norm(newray.direction)-1)<1e-8)
        else:
            #mirror
            newray.direction = ray.direction - 2*np.dot(ray.direction,normal)*normal        
        new_ray_bundle.append(newray)
    return new_ray_bundle

def partialsum(vector,ind):
    '''A function for partial sum'''
    return np.sum(vector[:(ind+1)])