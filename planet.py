from __future__ import division
from numpy import*
import random as r
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from scipy.constants import k
from AST1100SolarSystem import AST1100SolarSystem
import time

start = time.time()
system = AST1100SolarSystem(23078) #my seed

Mstar = system.starMass; Mplanets = system.mass #Mass of my plantes and star
init_x = system.x0; init_y = system.y0 #The initial position of my planets
init_vx = system.vx0; init_vy = system.vy0 # the initial velocity to my planets
Rplanets = system.radius; Rstar = system.starRadius #The radius of my plantes
"""
Comment to the last call (the radius), in case if I need it.
"""
g = 4*pi**2 #Newtons gravitational constant in astronomical units


class planet:

    def __init__(self,Mstar,Mplanets, Rplanets, Rstar, init_x,
    init_y, init_vx, init_vy,g):
        self.T = 30; self.time = 20000
        self.N = self.T*self.time;

        self.Mstar = Mstar; self.Mplanets = Mplanets
        self.Rplanets = Rplanets; self.Rstar = Rstar
        self.init_x = init_x; self.init_y = init_y
        self.init_vx = init_vx; self.init_vy = init_vy
        self.g = g; self.Nplanets = system.numberOfPlanets

        self.dt = 1/self.time


    def planet_position(self):
        vel =     zeros(shape=(2,self.Nplanets))
        pos =     zeros(shape=(2,self.Nplanets,self.N))
        a_p =     zeros(shape=(2,self.Nplanets))


        pos[0,:,0]  = self.init_x; pos[1,:,0] = self.init_y
        vel[0,:]    = self.init_vx; vel[1,:] = self.init_vy
        total_time = 0
        vel += 0.5*a_p*self.dt

        for i in xrange(0,self.N-1):
            for j in xrange(0,self.Nplanets):
                R = pos[:,j,i]
                a = (linalg.norm(R))**2
                b = (R/(linalg.norm(R)))

                a_p[:,j] = -((self.g*self.Mstar)/(a))*b

                vel[:,j] = vel[:,j] + a_p[:,j]*self.dt
                pos[:,j,i+1] = pos[:,j,i] + vel[:,j]*self.dt

        #system.checkPlanetPositions(pos,self.T,self.time)
        tid = linspace(0,20,20000*20)
        system.orbitXml(pos, tid)
        grid('on')
        xlabel('AU in dimension x')
        ylabel('Au in dimension y')
        title('Stationary star and the planets orbit')
        for m in xrange(0,self.Nplanets):
            plot(pos[0,m,:], pos[1,m,:])
            plot(pos[0,m,-1], pos[1,m,-1], '-o')
        plot(0,0, '*')
        #savefig('orbit_no_mov.png')
        show()

        #return pos, vel

    def big_planets(self, n=4):
        #index for the biggest planets
        bi = argpartition(self.Mplanets,-n)[-n:]
        return bi

    def star_position(self):
        vel =     zeros(shape=(2,self.Nplanets))
        pos =     zeros(shape=(self.N,2,self.Nplanets,))
        a_p =     zeros(shape=(2,self.Nplanets))


        pos[0,0,:]  = self.init_x; pos[0,1,:] = self.init_y
        vel[0,:]    = self.init_vx; vel[1,:] = self.init_vy
        #print pos

        a_star = zeros(shape = (2))

        a_pos = zeros(shape =(self.N,2))
        a_vel = zeros(shape =(self.N,2))

        total_time = linspace(0, self.N,self.N)

        for b in self.big_planets():
            a_vel[0,:] += -(self.Mplanets[b]*vel[:,b])/self.Mstar


        vel += 0.5*a_p*self.dt
        a_vel += 0.5*a_star*self.dt

        for i in xrange(0,self.N-1):
            for j in self.big_planets():
                R = pos[i,:,j] - a_pos[i,:]
                a = (linalg.norm(R))**2
                b = (R/(linalg.norm(R)))

                a_p[:,j] = -((self.g*self.Mstar)/(a))*b
                vel[:,j] = vel[:,j] + a_p[:,j]*self.dt
                pos[i+1,:,j] = pos[i,:,j] + vel[:,j]*self.dt

                ###############################
                Rr = pos[i,:,j] - a_pos[i,:]
                ar = (linalg.norm(Rr))**2
                br = (Rr/(linalg.norm(Rr)))

                a_star += ((self.g*self.Mplanets[j])/(ar))*br

            a_vel[i+1,:] = a_vel[i,:] + a_star*self.dt
            a_pos[i+1,:] = a_pos[i,:] + a_vel[i,:]*self.dt

            a_star = zeros(shape = (2))

              #for q in range(system.numberOfPlanets):
        grid('on')
        xlabel('AU in dimension x')
        ylabel('Au in dimension y')
        title('The stars movement when affected by big planets')
        for bigp in self.big_planets():
            plot(pos[:,0,bigp], pos[:,1,bigp])
            plot(pos[-1,0,bigp], pos[-1,1,bigp],'-o')

        plot(a_pos[:,0], a_pos[:,1])
        plot(a_pos[-1,0], a_pos[-1,1], '*')
        #savefig('planet_with_movstar.png')
        show()

        grid('on')
        xlabel('AU in dimension x')
        ylabel('Au in dimension y')
        title('The stars movement when affected by the big planets in the system')
        plot(a_pos[:,0], a_pos[:,1])
        plot(a_pos[-1,0], a_pos[-1,1], '*')
        #savefig('planet_with_movstar.png')
        show()

        grid('on')
        xlabel('time')
        title('velocity curve')
        ylabel('star velocity (AU)')
        plot(total_time, a_vel[:,0], '-r')
        #savefig('star_velocity_curve.png')
        show()

        print pos
        #system.checkPlanetPositions(pos,self.T,self.time)

        #return pos, vel

p = planet(Mstar,Mplanets, Rplanets, Rstar, init_x,
init_y, init_vx, init_vy,g)
p.planet_position()
#p.star_position()
end = time.time()
print ((end-start)), 'seconds'






####
