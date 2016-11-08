from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from AST1100SolarSystem import AST1100SolarSystem
import scipy.constants as scons
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
#from engine import motor
#from planet import planet
#from engine import engine_boost
import scipy.interpolate as inter
import time
start = time.time()

system = AST1100SolarSystem(23078) #my seed
#inFile = open('positionsHomePlanet.npy', 'rb')
planetPositions, times = np.load('planetPositions.npy')

class trajectory:

    def __init__(self, planetPositions, times):
        self.Mstar = system.starMass
        self.Mplanets = system.mass

        self.planetPositions = planetPositions
        self.init_x = system.x0
        self.init_y = system.y0

        self.init_vx = system.vx0
        self.init_vy = system.vy0

        self.home_m = system.mass[0]
        self.home_r = system.radius[0]

        self.Rplanets = system.radius
        self.Rstar = system.starRadius
        self.AU = 149597871*10**3

        self.Rau = self.home_r*10**3

        self.R = (self.Rau/self.AU)

        self.G = 4*np.pi**2
        self.time_step = 250000

        self.step = 10000000
        self.T = 10

        self.times = times


        self.N = self.T*self.time_step

        self.dt1 = 1/self.step
        self.dt2 = 1/self.time_step


        self.Nplanets = system.numberOfPlanets

        self.v_esc = np.sqrt((2*self.home_m*self.G)/(self.R))

        self.celestial_pos = inter.interp1d(self.times, self.planetPositions)
        #updating the celestial position ->star and planets

        self.star_pos = np.zeros(shape = (2,self.planetPositions[0,0,:].size))
        self.msat = 1100

    def masses(self): #A array containing the masses of the celestial objects
        mass = []
        for m in xrange(self.Nplanets):
            mass.append(self.Mplanets[m])
        mass.append(self.Mstar)
        all_mass = np.array(mass)

        return all_mass



    def deltav2(self,r1,r2): #to get into orbit at target planet
        mu = self.G*self.Mstar
        return np.sqrt(mu/r2)*(1-np.sqrt((2*r1)/r1+r2))


    def deltav1(self,r1,r2): #get out of orbit (homeplanet)
        mu = self.G*self.Mstar
        return np.sqrt(mu/r1)*(np.sqrt((2*r2)/(r1+r2))-1)*0.14838


    def shoot_angle(self,r1,r2): #The best angle to shoot out the rocket
        c = 1/(2*np.sqrt(2))
        f = (r1/r2 + 1)
        v = np.pi*(1-c*np.sqrt(f**3))

        return 180*v/np.pi


    def angle_bet_plan(self): # angle between the planetes, return time scalar
        #T = self.sat_period(self.semi_a(system.a[0],system.a[1]))
        #p = self.planet_period(system.a[1])

        theta = np.round(self.shoot_angle(system.a[0],system.a[1]))
        print theta
        angle_degrees = 0
        c = 0
        print "Please wait, so that the angle between targetplanet"
        print "and home planet is", theta, "degrees"
        print"----------------------------------------------------"

        while c < self.planetPositions[0,0,:].size:
            angle = np.arccos(np.dot(self.planetPositions[:,0,c],
            self.planetPositions[:,1,c])/(self.distance(0,c)*self.distance(1,c)))
            angle_degrees = int(180*angle/np.pi)
            t = self.times[c]
            c += 1
            if angle_degrees == theta:
                print "After",(t), "years, the angle between"
                print "homeplanet and targetplanet is", angle_degrees, "degrees"
                print "We are ready to launch"
                break
        print"----------------------------------------------------"

        launchtimer = 5
        print "This is spacecraft Poseidon, we are launching the rocket in"

        for i in range(5):
            launchtimer -= 1
            print launchtimer

        print "Rocket launched!"
        print "Poseidon this is Nasa, we wish you a nice trip!"

        return (t)


    def semi_a(self,r1,r2): #calculating the semi major axis between a p0 and p2
        return (r1 + r2)/2

    def sat_period(self,a): #the period it takes from p0 to p2
        return np.sqrt(a**3)/2

    def planet_period(self,a):
        return np.sqrt(a**3)


    def dist_plan(self,m):
        R = self.planetPositions[:,1,m] - self.planetPositions[:,0,m]
        return np.linalg.norm(R)

    def distance(self,n,m):
        return np.linalg.norm(self.planetPositions[:,n,m])

    def best_distance(self,k):
        return system.a[1]*np.sqrt((self.masses()[1]/self.Mstar)*1/k)

    def planet_velocity(self,planet, time):
        dtt = 10e-7

        return ((self.celestial_pos(time+dtt)[:,planet]-
        self.celestial_pos(time-dtt)[:,planet])/(2*dtt))

    def AUperYear_to_MeterPerSec(self):
        year = 60*60*24*365
        return self.AU/year

    def injection(self,time,planetnumber,rvec,vel_sat):
        mass = self.masses()
        theta = np.arctan2(rvec[1],rvec[0])
        orb_vel = np.sqrt(self.G*mass[1]/np.linalg.norm(rvec))

        return (np.array([orb_vel*np.sin(theta),-orb_vel*np.cos(theta)])-
        vel_sat + self.planet_velocity(1,time))



    def engine_info(self):
        mass_esc = 2.13e-13
        force =  1.75e-09

        return force, mass_esc
    def correction_burn(self,pnum,t):
        factor = self.planet_velocity(pnum,t)/np.linalg.norm(self.planet_velocity(pnum,t))
        return abs*factor

    def fuel(self,dv):
        unit = self.AUperYear_to_MeterPerSec()
        force_per_engine, numesc = self.engine_info()
        v = np.linalg.norm(dv) + 2.77
        return self.msat*(np.exp(((numesc*v)/force_per_engine)*unit)-1)

    def satelite_trajectory(self):
        launch_time= self.angle_bet_plan()
        pos_sat = np.zeros(shape = (2, self.N)) #satelite pos
        vel_sat = np.zeros(shape = (2)) #satelite vel

        q=0
        dv1 = self.deltav1(system.a[0],system.a[1])
        vel_sat= self.planet_velocity(0,launch_time)

        print (self.v_esc + dv1)*vel_sat/np.linalg.norm(vel_sat)


        vfac = self.v_esc + dv1
        vfac /= np.linalg.norm(vel_sat)
        vel_sat *= (1 + vfac)

        print 'initial velocity', vel_sat[0], vel_sat[1]


        pos_sat[0,q] = self.celestial_pos(launch_time)[0,0]
        pos_sat[1,q] = self.celestial_pos(launch_time)[1,0]
        pos_sat[:,q] += (self.R)*(vel_sat/np.linalg.norm(vel_sat))


        mass = self.masses()
        a_sat     = np.zeros(shape = (2))
        a_new_sat = np.zeros(shape = (2))
        print 'initialposition before takeoff is', pos_sat[0,q], pos_sat[1,q]


        #vel_sat += 0.5*a_sat*self.dt1
        time = launch_time
        list =[]
        ######## This section is launch section and out of home orbit###
        while time < launch_time + 0.06:
            a_sat[:] = 0
            R1 = -self.celestial_pos(time)[:,0] + pos_sat[:,q]
            R2 = pos_sat[:,q]

            a = (np.linalg.norm(R1))**3
            aa = (np.linalg.norm(R2))**3

            a_sat += -(((self.G*mass[0])/(a))*R1)
            a_sat = a_sat - ((self.G*self.Mstar)/aa)*R2
            vel_sat = vel_sat + a_sat*self.dt1
            pos_sat[:,q+1] = pos_sat[:,q] + vel_sat*self.dt1
            q += 1

            time += self.dt1
            list.append(time)
        ##################################################################

        print "-----------------------------------------------------------------"
        print " We have now left our homeplanets orbit, we are now headed to g23"
        print "-----------------------------------------------------------------"

        #vel_sat -= 0.5*a_sat*self.dt1
        #vel_sat += 0.5*a_sat*self.dt2

        t = time
        z = q

        rel_dist = []

        correct_burn = False
        burnedinorbit = False
        #time_correct = 10.75
        #### The probe is into deep space and is headed to planet 1#######
        for i in xrange(q,self.N-1):
            a_sat[:] = 0
            cel_pos = self.celestial_pos(t)
            for j in xrange(len(mass)-1):
                R1 = -cel_pos[:,j] + pos_sat[:,i]
                a = (np.linalg.norm(R1))**3
                a_sat += -(((self.G*mass[j])/(a))*R1)


            R2 = pos_sat[:,i]
            aa = (np.linalg.norm(R2))**3
            a_sat = a_sat - ((self.G*self.Mstar)/aa)*R2

            vel_sat = vel_sat + a_sat*self.dt2
            pos_sat[:,i+1] = pos_sat[:,i] + vel_sat*self.dt2

            sat_targetplan = pos_sat[:,i] - cel_pos[:,1]
            t += self.dt2
            #print np.linalg.norm(sat_targetplan)
            """if t < time_correct:
                if correct_burn == False:
                    print 'Our old velocity is', vel_sat, np.linalg.norm(vel_sat)
                    pos_sat[:,i] = np.array([-4.46381215281 ,  -5.20296058081])
                    vel_sat = np.array([6.43832690114 ,  -4.1614165334])
                    print 'Doing a correction burn, and readjusting velocity'
                    print 'New velocity', np.linalg.norm(vel_sat), vel_sat
                    correct_burn = True"""

            if np.linalg.norm(sat_targetplan) < 0.00047:
                #print np.linalg.norm(sat_targetplan)
                if burnedinorbit == False:
                    print 'old vel', np.linalg.norm(vel_sat),vel_sat
                    print 'old pos', pos_sat[:,i]
                    vel_sat += self.injection(t,1,sat_targetplan,vel_sat)
                    print 'injection boost', self.injection(t,1,sat_targetplan,vel_sat)
                    fuel_used = self.fuel(vel_sat)
                    print fuel_used, np.linalg.norm(vel_sat)
                    burnedinorbit = True
                    print 'doing a injection boost, to get in a stable orbit'
                    print 'Our new velocity is', np.linalg.norm(vel_sat)
                    print pos_sat[:,i], t, self.celestial_pos(t)[:,1]

            if  np.linalg.norm(sat_targetplan) < self.best_distance(10):
                rel_dist.append(sat_targetplan)



        ####################################################################
        return pos_sat,rel_dist

    def plot_trajectory(self):
        print 'Poseidon here, we are in orbit around planet g23'
        pos_sat, rel_dist = self.satelite_trajectory()

        for p in range(self.Nplanets):
            plt.plot(self.planetPositions[0,p,:], self.planetPositions[1,p,:])
            plt.plot(self.celestial_pos(5.98)[0,p], self.celestial_pos(5.98)[1,p], 'rh')
            plt.plot(self.planetPositions[0,p,0], self.planetPositions[1,p,0],'x')


        plt.grid('on')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Whole solarsystem with satellite orbit')
        plt.plot(pos_sat[0,:], pos_sat[1,:],'-k')
        plt.plot(pos_sat[0,0], pos_sat[1,0], '-kh')
        plt.plot(pos_sat[0,-1],pos_sat[1,-1], '-kh')
        plt.plot(0,0, '*')
        plt.savefig('whole_orbit.png')
        end = time.time()
        print (end-start), 'seconds!'
        plt.show()

        rel_dist = np.array(rel_dist)
        plt.grid('on')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('satellite orbit around G23')
        plt.plot(rel_dist[:,0], rel_dist[:,1])
        plt.plot(0,0, '-ro')
        plt.legend(['satellite (Poseidon)', 'G23'])
        plt.axis('equal')

        print (time.time()-start)
        plt.savefig('orbit_around_g23.png')
        plt.show()

    def plot_orbit(self):
        sat_pos, rel_dist = self.satelite_trajectory()

        rel_dist = np.array(rel_dist)
        plt.grid('on')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('satellite orbit around G23')
        plt.plot(rel_dist[:,0], rel_dist[:,1])
        plt.plot(0,0, '-ro')
        plt.legend(['satellite (Poseidon)', 'G23'])
        plt.axis('equal')

        print (time.time()-start)
        plt.savefig('orbit_around_g23.png')
        plt.show()

trajec = trajectory(planetPositions, times)
trajec.plot_trajectory()
#print trajec.injection(8.60840200009,1,np.array([-6.18241379,2.98485603]),np.array([2.36987566,-2.88066334]))
#trajec.angle_bet_plan()



####################
