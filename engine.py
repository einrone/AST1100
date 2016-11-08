from __future__ import division
from numpy import *
import random as r
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from scipy.constants import k
from AST1100SolarSystem import AST1100SolarSystem

myStarSystem = AST1100SolarSystem(23078)
massOfHomePlanet = myStarSystem.mass[0]
radiusOfHomePlanet = myStarSystem.radius[0]

##################################################################
Npart       = 100000
n           = 1000
L           = 1e-6    # nm, the length of the box
T           = 10000  # K, tempreture
total_time  = 1e-9
m           = 3.34*1e-27
sigma       = sqrt(k*T/m) #the mean velocity of the particles
dt          = total_time/(n-1)
home_m = massOfHomePlanet
home_r = radiusOfHomePlanet
#################################################################
class motor:

    def __init__(self,Npart,n,L,T,m,sigma,total_time,dt,home_m,home_r):
        self.Npart = Npart; self.n = n
        self.L = L; self.T = T; self.m = m
        self.sigma = sigma; self.dt = dt
        self.total_time = total_time
        self.home_r = home_r; self.home_m = home_m

    def box_collision_info(self):
        """
        In this function have given the particles the attribute to collide and
        and escape a hole within a interval L/4 <xy < 3L/4 when z < 0
        """
        position = np.zeros((self.Npart,3)) # antall part, dim, iterasjoner
        position[:,:] = np.random.uniform(0,1e-6, size = (self.Npart,3))
        velocity = np.zeros((self.Npart,3))
        velocity[:,:] = np.random.normal(0,self.sigma,size = (self.Npart,3))

        part_collided = 0
        part_escaped = 0
        momentum = 0

        print 'engine started'
        for i in xrange(1,self.n):
            #collision
            position += velocity*dt
            l_hole = position[:,0:2] > self.L/4
            h_hole = position[:,0:2] < (3*self.L)/4
            pos_xy = np.logical_and(l_hole, h_hole)
            pos_xy = np.logical_and(pos_xy[:,0], pos_xy[:,1])
            pos_z = position[:,2] < 0
            esc_part = np.logical_and(pos_z, pos_xy)

            #velocity[esc_part] = velocity[esc_part]
            part_escaped += np.sum(esc_part)

            for j in xrange(0,3):
                impact_wall_pos = np.logical_and(position[:,j] > 0,
                 position[:,j] < self.L)
                velocity[np.logical_not(impact_wall_pos),j] = -velocity[
                np.logical_not(impact_wall_pos),j]


                if j == 0:
                    part_collided += np.sum(np.logical_not(impact_wall_pos),j)
                    momentum += np.sum(2*self.m*abs(velocity[np.logical_not(
                    impact_wall_pos),j]))



            position[position < 0] = 0
            position[position >self.L] = self.L

        particle_collided = part_collided/2
        return position, velocity,part_escaped, impact_wall_pos, particle_collided, momentum



    def escaped_momentum(self):
        """
        calculates the escaped momentum per engine, because of this escape we
        get that the rocket accelerates upwards

        """
        position, velocity,escaped_particles,impact,collision,mom = self.box_collision_info()

        for i in xrange(1,self.n):
            velocity[np.logical_not(impact)] = velocity[np.logical_not(
            impact)]
            momentum = self.m*velocity
            abs_momentum = np.sum(np.sqrt(momentum[:,0]**2 + momentum[:,1]**2
            + momentum[:,2]**2))/2
        force = abs_momentum/self.dt

        return abs_momentum, force
    #def accelerartion_and_force(self):
    def pressure(self):
        momentum = 0
        pos, vel, part_esc, impact, collision, mom = self.box_collision_info()

        force = mom/(self.dt*(self.n-1))
        pressure = (force/self.L**2)/2
        return pressure


    def kinetic_energy(self):
        """
        calculates the invidual kinetic energy (by that i mean kinetic energy
        per particle), and it also calculate the total kinetic energy.
        """
        position, velocity, escaped_particles,impact, wall_collision,mom = self.box_collision_info()
        for j in xrange(1,self.n):
            abs_velocity = np.sqrt(velocity[:,0]**2+velocity[:,1]**2
            + velocity[:,2]**2)
            KE = 0.5*self.m*abs_velocity**2
            total_KE = np.sum(KE)
        invid_KE = total_KE/self.Npart

        return total_KE, invid_KE


    def checkFuelAmount(self, numberOfBoxes,fuel):
        esc_mom, per_force = self.escaped_momentum()
        pos, vel,part_esc, impact, part_coll, mom = self.box_collision_info()
        tot_force = self.engine_boost()
        gamma = 6.67e-11
        num_engine = 4.905e12
        esc_velocity = np.sqrt((2*gamma*self.home_m*2.0e30)/(self.home_r*10**3))


        seed = 23078
        myStarSystem = AST1100SolarSystem(seed)
        q = myStarSystem.massNeededCheck(numberOfBoxes, esc_velocity,
        tot_force, part_esc/(self.dt), fuel)



    def engine_boost(self):
        """
        This function calculates the force and acceleration per engine, by
        using the escaped momentum.

        """
        gamma = 6.67e-11
        esc_velocity = np.sqrt((2*gamma*self.home_m*2.0e30)/(self.home_r*10**3))
        print esc_velocity
        esc_mom, per_force = self.escaped_momentum()
        pos, vel,part_esc, impact, part_coll, mom = self.box_collision_info()

        rock_vel = []; rock_vel.append(0)
        rock_pos = []; rock_pos.append(0)
        #fuel = []; fuel.append(0)

        ##################constants####################
        ###############################################
        rocket_mass = 1100; num_engine = 1.186e13
        rocket_time = 20*60; up_time = 0
        mass_esc = (num_engine*part_esc*self.m)/self.total_time
        ###############################################

        delta_time = rocket_time/(1000)
        i = 0
        uptime = 0
        fuel = 55000
        velocity = 32100
        total_force = (esc_mom/self.total_time)*num_engine
        print total_force, 'newton'
        while (rock_vel[-1] < velocity and up_time < rocket_time and fuel > 0):
            total_acceleration = total_force/(rocket_mass + fuel)

            rock_vel.append(rock_vel[-1] + total_acceleration*delta_time)

            #fuel.append(fuel[-1] + mass_esc*delta_time)
            fuel -= mass_esc*delta_time

            i +=1
            uptime += delta_time


            if rock_vel[-1] > velocity:
                print "you have reached escape velocity"
                print rock_vel[-1]
                print fuel
                """myStarSystem.massNeededCheck(num_engine, esc_velocity,
                total_force/num_engine, part_esc/(self.total_time), fuel[-1])
                """
            if fuel < 0:
                print 'fuel done', fuel
                print rock_vel[-1], 'm/s'
                break

        plot(linspace(0,20*60, len(rock_vel)), rock_vel)
        show()
        return total_force


    def engine_and_general_info(self):
        """
        This function prints out important information about the engine.
        We se both analytical and numerical values.
        """
        pos,vel,esc_part, impact, wall_collision,mom = self.box_collision_info()
        tot_kin, kin_er = self.kinetic_energy()
        esc_mom, force = self.escaped_momentum()
        pres = self.pressure()
        tot_force = self.engine_boost()
        #force, acceleration, fuel = self.engine_boost()

        print"          Engine started and launched           "

        print "###############################################"
        print "      Engine status (Numerical values)         "
        print "-----------------------------------------------"
        print "The amount of particle escaped %g" %(esc_part)
        print "Amount of particles collided with one wall %i" %wall_collision
        print "Momentum escaped %g kgm/s" %(esc_mom)
        print "Kinetic energy per particle %gj" %(kin_er)
        print "Total kinetic energy %gj" %(tot_kin)
        print "Pressure inside the engine is %f" %(pres)
        print "momentum on the wall %g" %(mom)
        print "total force %g"%(tot_force)
        print "###############################################"
        print "             Launch info                       "
        print "-----------------------------------------------"
        #print "acceleration per engine %g m/s^2" %(acceleration)
        #print "force per engine %g N           " %(force)
        print "################################################"





k = motor(Npart,n,L,T,m,sigma,total_time,dt,home_m,home_r)
k.engine_and_general_info()
















#
