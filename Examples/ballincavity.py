import matplotlib.pyplot as plt
from pylab import *
from scipy.optimize import brentq
import numpy as np


class Ray:
    def __init__(self):
        self.direction = array([0,0,1])
        self.r0 = array([0,0,0])
        self.length = 0
        self.wavelength=589
    def __str__(self):
        return "Ray starting at "+str(self.r0)
    def plot(self):
        r0 = self.r0
        l = self.length;
        r1 = (self.r0 + l*self.direction )
        plt.plot((r0[2],r1[2]),(r0[0],r1[0]),'b')



class Ball:
    #the inside is defined as the region F>0
    def __init__(self):
        self.center = array([0,0,0])
        self.radius = 1
        self.n_inside = 1.0
        self.n_outside = 1.8
        self.label="round thing"
    def __str__(self):
        return "Ball with label %s"%(self.label)

    def plot(self):
        thetavals = np.linspace(0,2*np.pi,1024)
        zvals = self.center[2]+self.radius*np.sin(thetavals)
        xvals = self.center[0]+self.radius*np.cos(thetavals)
        plt.plot(zvals,xvals,'r')
    
    def Ffunc(self,r):
        #returns radius^2 - |r-center|^2
        return self.radius**2 - dot(r-self.center,r-self.center)
    
    def gradFfunc(self,r):
        return 2*(self.center-r)
                
    def find_intersect(self,ray): # does the ray intercept the surface if so how far along ray? if not return -1
        s0 = 1e-5 #start this distance along the ray so that we don't get confused by the surface we have just passed thru
        
#        print(ray.r0)
#        print(self.Ffunc(ray.r0))
        if self.Ffunc(ray.r0+s0*ray.direction)>0: #do we start of inside the sphere
            #code for if we do:
            #if you go three radii in any direction you are definatley outside the sphere so this point is definately outside the sphere
            soutside = 3*self.radius
            #now we have values of s that bracket the intercept we can find it
            s = brentq(lambda x: self.Ffunc(ray.r0+x*ray.direction),s0,soutside)
            return s
        else:
            #code for if we start outside sphere
            #s for closest point to the center
            s_closest = dot(self.center-ray.r0,ray.direction)
            if s_closest < s0:
                return -1 #we are outside and heading away from center
            #f value for closest point to the center
            F_closest = self.Ffunc(ray.r0+s_closest*ray.direction)
            if F_closest<0: #ray never hits sphere
                return -1
            else:
                #s_closest and 0 are on different sides of the sphere
                #print(s_closest)
                #print(ray.r0)
                #print(ray.direction)
                #print(F_closest)
                #print(self.Ffunc(ray.r0))
                s = brentq(lambda x: self.Ffunc(ray.r0+x*ray.direction),0,s_closest)
                return s

        #we shouldn't never reach this point
        #this is a way to make sure know if we do
        assert(False) 



        
def trace(ray_bundle,surfs):
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
#        print("smin=",smin)
            

        ray.length=smin
        r = ray.r0+(smin)*ray.direction
        normal = nextsurf.gradFfunc(r)
        normal = normal/norm(normal)

        newray = Ray()
        newray.r0 = r
        newray.length=0.2


        #this chooses the direction of the new ray
        n1 = nextsurf.n_outside
        n2 = nextsurf.n_inside
        # if either of the refractive indicies are infinite
        # then act like is't s a mirror
        # otherwise use snells law
        if np.isfinite(nextsurf.n_outside) and np.isfinite(nextsurf.n_inside):
            # donesn't look like snells law
            # but see the a paper
            if dot(ray.direction,normal)<0:
                n2 = nextsurf.n_outside
                n1 = nextsurf.n_inside
            else:
                n1 = nextsurf.n_outside
                n2 = nextsurf.n_inside
                
            a = n1/n2*dot(ray.direction,normal)
            b = ((n1/n2)**2-1)
            gamma = -b/(2*a)
            for k in range(30):
                #        print k,gamma
                gammanew = (gamma**2-b)/(2*(gamma+a))
                if abs(gamma-gammanew) <1e-8:
                    gamma = gammanew
                    break
                gamma=gammanew
            #make sure iterative approach has converged 
            assert(k<25)
            newray.direction = n1/n2*ray.direction+gamma*normal
            assert(abs(norm(newray.direction)-1)<1e-8)
        else:
            #mirror
            newray.direction = ray.direction - 2*dot(ray.direction,normal)*normal        
        new_ray_bundle.append(newray)



        
    return (new_ray_bundle)
   
#################################################
## setup the problem
#######################################

plt.close('all')
plt.ion()

b = Ball()
b.n_inside=1.6
b.n_outside=1
b.label = "centralball"

m1 = Ball()
m1.radius = 20
m1.n_inside=1
m1.n_outside=np.inf #it's a mirror
m1.center = array([0,0,10])
m1.label="mirror1"

m2 = Ball()
m2.radius = 20
m2.center = array([0,0,-10])
m2.n_inside=1
m2.n_outside=np.inf #it's a mirror
m2.label="mirror2"
b.plot()
m1.plot()
m2.plot()




rays = []

#for x in linspace(.75,-.75,5):
#     r = Ray()
#     r.r0 = array([x,0,-2])
#     r.length = 1
#     rays.append(r)

for th in linspace(-0.1,0.1,5):
    r = Ray()
    r.r0 = array([0,0,-3])
    r.direction = array([np.sin(th),0,np.cos(th)])
    r.length = 1
    rays.append(r)
    


surfs = [b,m1,m2]
for k in range(4):
    newrays = trace(rays,surfs)
    for r in rays:
        r.plot()
    rays=newrays
    
         
     
     
plt.axis((-11,11,-3,3))
#plt.axis('image')
plt.show()
     
