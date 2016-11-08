from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from AST1100SolarSystem import AST1100SolarSystem
import scipy.constants as scons
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import scipy.interpolate as inter
from PIL import Image
import sys
import time
start = time.time()

system = AST1100SolarSystem(23078) #my seed
planetPositions, times = np.load('planetPositions.npy')



class skynet:

    def __init__(self):
        self.ypix = 480
        self.xpix = 640
        self.FOWtheta = (np.pi*70)/180 #field of wiew in theta
        self.FOWphi = (np.pi*70)/180 #field of wiew in phi


        inFile = open('himmelkule.npy', 'rb')
        self.himmelkulen = np.load(inFile)
        inFile .close()

        self.times = times
        self.planetPositions = planetPositions
        self.Nplan = system.numberOfPlanets
        self.intplanpos = inter.interp1d(self.times,self.planetPositions)

    def max_min(self):
        xmin = - (2*np.sin(self.FOWphi/2))/(1+np.cos(self.FOWphi/2))
        xmax =   (2*np.sin(self.FOWphi/2))/(1+np.cos(self.FOWphi/2))

        ymin = - (2*np.sin(self.FOWtheta/2))/(1+np.cos(self.FOWtheta/2))
        ymax =   (2*np.sin(self.FOWtheta/2))/(1+np.cos(self.FOWtheta/2))

        return xmin, xmax, ymin, ymax

    def degrees_to_radians(self,angle):
        return (np.pi*angle)/180

    def radians_to_degrees(self,radians):
        return (180*radians)/np.pi

    def picture(self,theta0, phi0):
        xmin, xmax, ymin, ymax = self.max_min()

        x_picture = np.linspace(xmin, xmax, self.xpix)
        y_picture = np.linspace(ymax, ymin, self.ypix)

        picture = np.zeros((self.ypix, self.xpix,3),dtype = np.uint8)
        """
        for every time we calculate ypic, we have to calculate 640 points in
        xpic direction.
        """

        for i in xrange(self.xpix):
            #y picture direction, we go from 0 to 480
            for j in xrange(self.ypix):
            # x picture direction, we got from 0 640

                rho = np.linalg.norm(np.array([x_picture[i], y_picture[j]]))
                c = 2*np.arctan2(rho,2)

                eq1 = np.cos(c)*np.cos(theta0) + (
                (y_picture[j]*np.sin(c)*np.sin(theta0))/rho)

                theta = np.pi/2 - np.arcsin(eq1)

                eq2 = (x_picture[i]*np.sin(c))/(
                rho*np.sin(theta0)*np.cos(c)-y_picture[j]*np.cos(theta0)*np.sin(c))

                phi = phi0 + np.arctan(eq2)


                pixnum = AST1100SolarSystem.ang2pix(theta, phi)
                temp = self.himmelkulen[pixnum]
                rgb = [temp[2] , temp[3] , temp[4]]
                picture[j,i,:] = rgb

        return picture

    def picture_generator(self,name):
        picture = np.load('orient0.npy')#np.load('projection_image.npy')

        try:

            img2 = Image.fromarray(picture[190,:,:,:])
            img2.save(name + '.png')
        except ValueError:
            print ('You did not typa string, please enter a string')
            img2 = Image.fromarray(picture[0,:,:,:])
            img2.save(name + '.png')


    def projection(self,theta0,phi0,name):

        projection_image = np.zeros((360,self.ypix, self.xpix, 3), dtype = np.uint8)

        for i in xrange(360):
            phi = self.degrees_to_radians(i)
            projection_image[i,:,:,:] = self.picture(theta0,phi0)

        #outfile = projection_image
        np.save(name + '.npy', projection_image)


    def compare(self,projectionfile, satelitepicture):
        picture360 = projectionfile
        takenpic = Image.open(satelitepicture)
        angle = 360
        phi = 0
        total_error = 10000

        for i in xrange(angle):
            error = np.sum((takenpic - picture360[i,:,:,:])**2)
            if error < total_error:
                phi = i
                total_error = error
        return phi




    def values(self):
        dL1 = 0.003497832158; phi1 = self.degrees_to_radians(300.499476)
        dL2 = -0.003149620043; phi2 = self.degrees_to_radians(244.691555)
        L0 = 656.3

        return dL1,dL2, phi1, phi2, L0

    def velocityrefstar(self):
        dL1, dL2, phi1, phi2, L0 = self.values()
        c = self.speed_of_light_in_au()

        vsref1 = (dL1/L0)*c
        vsref2 = (dL2/L0)*c
        return vsref1, vsref2

    def speed_of_light_in_au(self):
        C_au= 63239.7263
        return C_au


    def velocitysat(self,dl1_rel, dl2_rel):
        dL1, dL2, phi1, phi2, L0 = self.values()
        c = self.speed_of_light_in_au()
        vsref1, vsref2 = self.velocityrefstar()

        vsrel1 = (dl1_rel/L0)*c
        vsrel2 = (dl2_rel/L0)*c

        v1 = vsref1 - vsrel1
        v2 = vsref2 - vsrel2


        Vx = (1/(np.sin(phi2-phi1)))*(np.sin(phi2)*v1 - np.sin(phi1)*v2)
        Vy = (1/(np.sin(phi2-phi1)))*(-np.cos(phi2)*v1 + np.cos(phi1)*v2)

        return Vx, Vy

    def sat_position(self,distance, time, i, j):

        x1,y1 = self.intplanpos(time)[:,i]
        x2,y2 = self.intplanpos(time)[:,j]

        diststar = distance[-1] #star
        dist1 = distance[i] #planet i
        dist2 = distance[j] #other planet j

        a = np.array(([x1,y1],[x2,y2]))
        b = np.array([(diststar**2 - dist1**2 + x1**2 + y1**2)*0.5,
                     (diststar**2 - dist2**2 + x2**2 + y2**2)*0.5])

        satelite_position = np.linalg.solve(a,b)

        return satelite_position



#print system.getRefStars()

projectionfile = np.load('projection_image.npy')
posfile = np.load('pos.npy')
#print posfile
p = skynet()
#print p.velocitysat(0.129494449336,0.070932854259)
p.picture_generator()
np.load('orient0.npy')
end = time.time()
print (end-start)
##########
