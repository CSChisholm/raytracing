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

        assert(false) #code should never reach here


    def bend(self,ray,s):
        #given a ray and the value of s where it intersects the surface
        #trim the ray and make a new bent one
        ray.length = s
        
def trace(ray_bundle,surfs):
    new_ray_bundle=[]

    for ray in ray_bundle:
        #work out the next surface for the ray to hit
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

        n1 = nextsurf.n_outside
        n2 = nextsurf.n_inside
        if np.isfinite(n1) and np.isfinite(n2):
            #snells law!
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

            assert(k<25)
            newray.direction = n1/n2*ray.direction+gamma*normal
            assert(abs(norm(newray.direction)-1)<1e-8)
        else:
            #mirror
            newray.direction = ray.direction - 2*dot(ray.direction,normal)*normal
                    
            

        new_ray_bundle.append(newray)



        
    return (new_ray_bundle)

#        normal = surf.gradFfunc(r)
#        normal = normal/norm(normal)
        
    
    #calculates the center of curvature of the surface
#     if abs(surf.curv)>1/100:
#         ccurv = array([0,0,surf.Z0+1./surf.curv])
#         print(ccurv)
#     else:
#         print('nope\n\n')
#         ccurv = None

#     count=0
#     for ray in ray_bundle:
#     #first determine value for s where the ray crosses the z= z0 plane
#         #ccurv for the lens
#         count=count+1
#         print("Ray = "+str(count))
#         if type(ccurv)!=type(None):
#             tv = (ray.r0-ccurv) #temp vector
#             b = 2*dot(ray.direction,tv)
#             c = dot(tv,tv)-1/surf.curv**2
#             s1 = (-b -sqrt(b**2-4*c))/2
#             s2 = (-b+sqrt(b**2-4*c))/2
#             if surf.curv >0:
#                 s = s1
#             else:
#                 s = s2
#             print(s1, s2, s,b**2-4*c)
#             if b**2-4*c<0:
#                 s = (ray.r0[2]-surf.Z0)/ray.direction[2]
#         else:
#             s = (ray.r0[2]-surf.Z0)/ray.direction[2]
#         for k in range(30):
#             r = ray.r0+s*ray.direction
#             dFds = dot(ray.direction,surf.gradFfunc(r))
#             snew = s - surf.Ffunc(r)/dFds
#             if abs(s-snew)<1e-9:
#                 break
#             s = snew
#         if k>25:
#             #assert(k<25)
#             #here we assume that the ray has missed the surface so for the new ray make r0=r0 and length of old ray zero
#             ray.length=0
#             newray = Ray()
#             newray.r0 = ray.r0
#             newray.direction = ray.direction
#         else:
#             ray.length=snew
#             r = ray.r0+snew*ray.direction
#             normal = surf.gradFfunc(r)
#             normal = normal/norm(normal)
        
#             newray = Ray()
#             newray.r0 = r
#             print(dot(ray.direction,normal))
#             a = surf.n1/surf.n2*dot(ray.direction,normal)
#             b = ((surf.n1/surf.n2)**2-1)
#             gamma = -b/(2*a)
#             for k in range(30):
#                 #        print k,gamma
#                 gammanew = (gamma**2-b)/(2*(gamma+a))
#                 if abs(gamma-gammanew) <1e-8:
#                     gamma = gammanew
#                     break
#                 gamma=gammanew
    
#             assert(k<25)
# #    print gamma,a,b
# #    print gamma**2+2*a*gamma+b
# #    print normal,norm(normal)
    
#             newray.direction = surf.n1/surf.n2*ray.direction+gamma*normal
        


#         assert(abs(norm(newray.direction)-1)<1e-8)
       
#         new_ray_bundle.append(newray)
    
#     return (new_ray_bundle)
    


plt.close('all')
plt.ion()

b = Ball()
b.n_inside=1
b.n_outside=1
b.label = "centralball"

m1 = Ball()
m1.radius = 20
m1.n_inside=1
m1.n_outside=np.inf
m1.center = array([0,0,15])
m1.label="mirror1"

m2 = Ball()
m2.radius = 20
m2.center = array([0,0,-15])
m2.radius = 20
m2.n_inside=1
m2.n_outside=np.inf
m2.label="mirror2"
b.plot()
m1.plot()
m2.plot()

rays = []
for x in [0.5]:#linspace(.75,-.75,5):
     r = Ray()
     r.r0 = array([x,0,-2])
     r.length = 1
     rays.append(r)


surfs = [b,m1,m2]
for k in range(500):
    newrays = trace(rays,surfs)
    for r in rays:
        r.plot()
    rays=newrays
    
         
     
     
plt.axis((-7,7,-3,3))
#plt.axis('image')
plt.show()
     



        
# class SphericalSurface:
#     def __init__(self):
#         self.Z0 = 0
#         self.curv = 0 
#         self.n1=1.
#         self.n2=2.
#         self.app=6.

#     def Ffunc(self,r):
#         ccurv = [0,0,self.Z0+1/self.curv]
#         F = dot(r-ccurv,r-ccurv)-(1/self.curv)**2
#         return F
#     def gradFfunc(self,r):
#         ccurv = [0,0,self.Z0+1/self.curv]
#         return 2*(r-ccurv)
#     def Zfunc(self,ssq):
#         #not used in calcs
#         c = self.curv
#         Z = self.Z0 + c*ssq/(1+sqrt(1-c*c*ssq))
#         return Z    


# class Surface:
#     def __init__(self):
#         self.Z0 = 0
#         self.curv = 0
#         self.kappa = 1 #=1+k
#         self.Aparams = array([])
#         self.n1 = 1.0
#         self.n2 = 1.0
#         self.app = 6
#     def Zfunc(self,ssq):
#         c = self.curv
#         Z = self.Z0 + c*ssq/(1+sqrt(1-self.kappa*c*c*ssq))
#         for k in range(len(self.Aparams)):
#             Z = Z+self.Aparams[k]*ssq**(2+k)
#         return Z    
            
#     def Ffunc(self,r):
#         x = r[0]
#         y = r[1]
#         z = r[2]
#         ssq = x*x+y*y
#         F = z-self.Zfunc(ssq)
#         return F

#     def gradFfunc(self,r):
#         x = r[0]
#         y = r[1]
#         z = r[2]
#         ssq = x*x+y*y
#         c = self.curv
#         E = c/sqrt(1-self.kappa*c*c*ssq)
#         for k in range(len(self.Aparams)):
#             E = E + 2*(2+k)*self.Aparams[k]*ssq**(1+k)#this could be wrong
#         F = z-self.Zfunc(ssq)
#         Fx = -x*E
#         Fy = -y*E
#         Fz = 1 
#         return array([Fx,Fy,Fz])
        


# def trace(ray_bundle,surf):
#     new_ray_bundle=[]

    
#     #calculates the center of curvature of the surface
#     if abs(surf.curv)>1/100:
#         ccurv = array([0,0,surf.Z0+1./surf.curv])
#         #print(ccurv)
#     else:
#         #print('nope\n\n')
#         ccurv = None

#     count=0
#     for ray in ray_bundle:
#     #first determine value for s where the ray crosses the z= z0 plane
#         #ccurv for the lens
#         count=count+1
#         #print("Ray = "+str(count))
#         if type(ccurv)!=type(None):
#             tv = (ray.r0-ccurv) #temp vector
#             b = 2*dot(ray.direction,tv)
#             c = dot(tv,tv)-1/surf.curv**2
#             s1 = (-b -sqrt(b**2-4*c))/2
#             s2 = (-b+sqrt(b**2-4*c))/2
#             if surf.curv >0:
#                 s = s1
#             else:
#                 s = s2
#             #print(s1, s2, s,b**2-4*c)
#             if b**2-4*c<0:
#                 s = (ray.r0[2]-surf.Z0)/ray.direction[2]
#         else:
#             s = (ray.r0[2]-surf.Z0)/ray.direction[2]
#         for k in range(30):
#             r = ray.r0+s*ray.direction
#             dFds = dot(ray.direction,surf.gradFfunc(r))
#             snew = s - surf.Ffunc(r)/dFds
#             if abs(s-snew)<1e-9:
#                 break
#             s = snew
#         if k>25:
#             #assert(k<25)
#             #here we assume that the ray has missed the surface so for the new ray make r0=r0 and length of old ray zero
#             ray.length=0
#             newray = Ray()
#             newray.r0 = ray.r0
#             newray.direction = ray.direction
#         else:
#             ray.length=snew
#             r = ray.r0+snew*ray.direction
#             normal = surf.gradFfunc(r)
#             normal = normal/norm(normal)
        
#             newray = Ray()
#             newray.r0 = r
#             #print(dot(ray.direction,normal))
#             a = surf.n1/surf.n2*dot(ray.direction,normal)
#             b = ((surf.n1/surf.n2)**2-1)
#             gamma = -b/(2*a)
#             for k in range(30):
#                 #        #print k,gamma
#                 gammanew = (gamma**2-b)/(2*(gamma+a))
#                 if abs(gamma-gammanew) <1e-8:
#                     gamma = gammanew
#                     break
#                 gamma=gammanew
    
#             assert(k<25)
# #    #print gamma,a,b
# #    #print gamma**2+2*a*gamma+b
# #    #print normal,norm(normal)
    
#             newray.direction = surf.n1/surf.n2*ray.direction+gamma*normal
        


#         assert(abs(norm(newray.direction)-1)<1e-8)
       
#         new_ray_bundle.append(newray)
    
#     return (new_ray_bundle)


    

# plt.close('all')
# rays = []


# for x in linspace(1.5,-1.5,20):
#     r = Ray()
#     r.r0 = array([x,0,-8.5])
#     rays.append(r)


# #
# #  set surfaces
# # 


# surfs = []



# #parameters of the lightpath 352330 aspheric
# as1curv = 1/-3.20561
# as1kappa = 1+-12.418013
# as1Aparams = array([9.00531e-3,-1.359752e-3,1.136638e-4,-4.278925e-6])
# athickness=3.19
# awd = 1.76
# as2curv = 1/2.74797
# as2kappa = 1+-0.542698
# as2Aparams = array([-3.19546e-4,-4.397785e-5,1.842256e-5,-1.566446e-6])


# nglass = 1.595
# nyso = 1.8
# ballradius = 0.625

# zshort = ballradius/nyso
# zlong = ballradius*nyso

# offset = +zlong-zshort


# s1 = Surface()
# s1.curv = as2curv
# s1.kappa = as2kappa
# s1.Aparams = as2Aparams
# s1.Z0 = -1.76-3.19+offset
# s1.n2 = nglass




# s2 = Surface()
# s2.curv = as1curv
# s2.kappa = as1kappa
# s2.Aparams = as1Aparams
# s2.Z0 = -1.76+offset
# s2.n1 = nglass


# s3 = SphericalSurface()
# s3.curv = 1/ballradius
# s3.Z0 = -ballradius-zshort
# s3.n2=nyso
# s3.app = 2*ballradius

# #s4 = SphericalSurface()
# #s4.curv = -1/ballradius
# #s4.Z0 = ballradius+zshort
# #s4.n1=nyso
# #s4.app = 2*ballradius

# #s5 = Surface()
# #s5.curv = -as1curv
# #s5.kappa = as1kappa
# #s5.Aparams = -as1Aparams
# #s5.Z0 = 1.76-offset
# #s5.n2 = nglass


# #s6 = Surface()
# #s6.curv = -as2curv
# #s6.kappa = as2kappa
# #s6.Aparams = -as2Aparams
# #s6.Z0 = 1.76+3.19-offset
# #s6.n1 = nglass

# s7=Surface()
# #s7.Z0=8.5
# s7.Z0=0.01

# surfs.append(s1)
# surfs.append(s2)
# surfs.append(s3)
# #surfs.append(s4)
# #surfs.append(s5)
# #surfs.append(s6)
# surfs.append(s7)



# for s in surfs:
#     xvals = linspace(-s.app,s.app,100)/2
#     zvals = 0*xvals;
#     for k in range(len(xvals)):
#         zvals[k] = s.Zfunc(xvals[k]**2)
#     plt.plot(zvals,xvals,'r')

#     newrays = trace(rays,s)
#     for r in rays:
#         r.plot()
       
#     rays=newrays

# for r in rays:
#     r.plot()

# plt.axis([-5.5,5.5,-3,3])    
# plt.axis('equal')
# plt.show()



    







        
        
    

                        
                        

                        

                        
                       
    
