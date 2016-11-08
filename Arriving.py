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

class lander:

    def __init__(self):

        self.T0 = 273.38; self.kb = 1.38065e-23
        self.mu = 30.114; self.mH = 1.67e-27; self.G = 6.67e-11
        self.gamma = 1.4

        self.a1 = system.a[1]; self.Rstar = system.starRadius
        self.a0 = system.a[0]
        self.Tstar = system.temperature; self.rho0 = system.rho0[1]

        self.Mp = system.mass[1]; self.Rp  = system.radius[1]
        self.num = ((self.gamma-1)/self.gamma)*self.mu*self.mH*self.Mp*self.G
        self.P0 = 0;

        self.denominator = self.kb*self.T0
        self.factor = self.num/self.denominator
        self.expo1 = self.gamma/(self.gamma-1)
        self.expo2 = 1/(1-self.gamma)

        self.omega = (self.kb*self.T0*self.Rp)/(2*self.mu*self.mH*self.Mp*self.G)
        self.R0 = (self.Rp*self.omega)/(1 + self.omega)
        self.Msat = 90 #kg
        self.Cd = 1
        self.N = 100000
        self.period = 90*self.N
        self.dt = 1/self.N
        self.yr = 60*60*24*365
        self.zz = self.G*self.Mp
        self.mult = np.log(1e-12/self.rho0)
        self.beta = (self.kb*self.T0/(self.mu*self.mH*self.Mp*self.G))
        self.d = -self.beta*self.Rp/(1+self.beta)

    def dv2(self):

        hohmann = np.sqrt(self.zz/self.a1)*(1-np.sqrt((2*self.a0)/(
        self.a0 + self.a1)))

        return hohmann




    def atmosphere_temperature(self,R):
        if z > self.R0: #isothermal atmosphere
            return self.T0/2

        if z <= self.R0: #adiabatic atmosphere
            return self.T0*(1-self.factor*(1/(R+self.Rp)-1/self.Rp))

    def atmosphere_pressure(self,R):
        if z <= self.R0: #adiabatic area
            return self.P0*(1-self.factor(1/(R+self.Rp) - 1/self.Rp))**(self.expo1)
        if z > self.R0: #isothermal area
            return self.P0*np.exp(-(2*self.mu*self.mH*self.Mp*self.G*R)/(
            self.kb*self.T0*(R + self.Rp)))

    def atmosphere_denisty(self,R):
        if z <= self.R0: #adiabatic area
            return self.rho0*(1-self.factor*(1/(R+self.Rp)-1/self.Rp))**(self.expo2)
        if z > self.R0: #isothermal area
            return  self.rho0*np.exp(-(2*self.mu*self.mH*self.Mp*self.G*R)/(
            self.kb*self.T0*(R + self.Rp)))

    def dragfriction(self,R,v):
        return 0.5*self.atmosphere_denisty(R)*self.Cd*self.A*v**2

    def parachute(self):
        vt = 3
        return (vt**2*self.Rp*self.rho0)/(2*self.G*self.Mp*self.Msat)

    def ratio(self,Fg,Fd):
        return Fg/Fd


    def terminal_velocity(self,A):
        return np.sqrt((2*self.G*self.Msat*self.Mp)/(self.Rp*A*self.rho0))

    def Fg(self,r,rvec):
        return ((self.G*self.Mp*self.Msat)/(np.linalg.norm(r))**3)*rvec

    def velocity_vec_0(self):
        vvec = np.array([-1.49009606085 , -3.65850904163,0 ])
        absvvec = np.linalg.norm(vvec)
        return absvvec

    def loworbit(self):
        pos_sat = np.zeros((3,self.period))
        vel_sat = np.zeros((3,self.period))

        velocity_vector = self.velocity_vec_0()
        hohmann = self.dv2()
        dv = velocity_vector*hohmann

        pos_sat[:,0] = np.array([ -6.3391849735 , 2.60854562701,0])
        vel_sat[:,0] = np.array([-1.49009606085 , -3.65850904163,0 ]) + dv

        t = 8.71*self.yr

        for i in xrange(self.period-1):
            rvec = pos_sat[:,i]
            a = self.Fg(pos_sat[:,i],rvec)
            vel_sat[:,i+1] = vel_sat[:,i] + a*self.dt
            pos_sat[:,i+1] = pos_sat[:,i] + vel_sat[:,i+1]*self.dt
            t += self.dt

            if np.linalg.norm(pos_sat[:,i]) < self.d:
                print 'we are near the G23 atmosphere'
                print 'time', t, 'at position', pos_sat
                plt.plot(pos_sat[0,:], pos_sat[1,:],'-')
                plt.plot(0,0, 'or')
                plt.show()
                break
            if t > 8.710003*self.yr:
                print 'bro you gotta do it better'
                plt.plot(pos_sat[0,:], pos_sat[1,:],'-')
                plt.plot(0,0, 'or')
                plt.show()
                break





p = lander()
p.loworbit()
