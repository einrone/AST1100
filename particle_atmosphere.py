from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from AST1100SolarSystem import AST1100SolarSystem
import scipy.constants as scons
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import scipy.interpolate as inter
import time
start = time.time()

system = AST1100SolarSystem(23078)
class atmosphericparticles:

    def __init__(self, Filespectrum, Filesigma):
        #constants and different variables/parameters
        self.spectrum = np.load(Filespectrum + '.npy')
        self.sigmanoise = np.load(Filesigma + '.npy')

        self.a = system.a[1]
        self.e = system.e[1]
        self.Ts = system.temperature
        self.Rstar = system.starRadius
        self.Ms = system.starMass
        self.AU = 149597871 #km
        self.Rs = self.Rstar/self.AU
        self.diviationN = self.sigmanoise[:,1]


    def files(self,spectrumfile, sigmafile): #makes the .txt file into .npy file
        wenttrough1 = False
        wenttrough2 = False

        if spectrumfile.endswith('.txt'):
            if wenttrough1 == False:
                spectratxt = np.loadtxt(spectrumfile)
                specfile = np.save('spectrumfile.npy', spectratxt)
                wenttrough1 = True

        if sigmafile.endswith('.txt'):
            if wenttrough2 == False:
                sigmatxt = np.loadtxt(spectrumfile)
                sigmafile = np.save('sigmafile.npy',sigmatxt)
                wenttrough2 = True


        else:
            print 'please type in a text file'

    def T0max(self,a,e,Ts,Rs): #maximum effective temperature, stefan-boltz law
        return Ts*np.sqrt(Rs/(2*a*(1-e**2)))

    def T0min(self,a,Ts,Rs):#minimum effective temperature, stefan-boltz law
        return Ts*np.sqrt(Rs/(2*a))

    def sentripetal_velocity(self): #the satelites velocity relative to G23
        v = 10000 #m/s
        return v

    def gases(self): # gases with their absorption lines

        O2 = np.array([630,690,760]) #oxygen absoption spectre 0
        CH4 = np.array([1660,2200]) #methane absoption spectre 1
        NO = np.array([2870]) #nitrous oxide absoption spectre 2
        CO = np.array([2340]) #carbon monooxide absoption spectre 3
        CO2 = np.array([1440,1600]) #carbon dioxide absoption spectre 4
        H2O = np.array([720,820,940]) #water vapour absoption spectre 5

        gas = np.array([O2,CH4,NO, CO, CO2, H2O])
        return gas

    def gasmass(self):
        unit = 1.660539040e-27
        mass = [31.9988, 16.04246, 30.0064, 28.0094,48.088,18.0144]
        mass = np.array(mass)*unit

        return mass

    def typegas(self):
        O2 = np.array(['O2_630nm', 'O2_690nm', 'O2_760nm'])
        CH4 = np.array(['CH4_1660','CH4_2200'])
        NO = np.array(['NO_2870nm'])
        CO = np.array(['CO_2340nm'])
        CO2 = np.array(['CO2_1440nm','CO2_1600nm'])
        H2O = np.array(['H2O_720nm','H2O_820nm','H2O_940nm'])

        return np.array([O2,CH4,NO, CO, CO2, H2O])

    def FWHM(self,lmbdcenter,T,m): #Full width at half maximum
        c = 299792458 #m/s
        k = 1.3806505e-23 #boltszmann constant
        return (2*lmbdcenter/c)*np.sqrt((2*k*T*np.log(2))/m)

    def width(self,FWHM): #Finds the width of the graph, sigma. This is used to
        #find F(lambda) in the next function
        return FWHM/(2*np.sqrt(2*np.log(2)))

    def flux(self,lmbd,Fmin,width, lmbdC): #gaussian flux
        Fmax = 1
        exponential = np.exp(-(lmbd - lmbdC)**2/(2*width**2))
        return Fmax + (Fmin - Fmax)*exponential

    def satlambda(self,lmbdC): # Observed wavelength
        velrel  = self.sentripetal_velocity() #satelites velocity

        ltow    = lmbdC + (velrel/(3e8))*lmbdC #positive direction
        lback   = lmbdC - (velrel/(3e8))*lmbdC #negative direction

        return ltow, lback
    def plot_lines(self,wl,obsflux,fluxmodel,a):
        #typegases = self.typegas()

        plt.plot(wl, obsflux, '-g')
        plt.hold('on')
        plt.plot(wl,fluxmodel,'-k')
        plt.title(str(a) + 'nm')
        plt.savefig(str(a) + 'nm' + '.png')
        plt.show()

    def interesting_lines_observed(self):
        typegases = self.typegas() #name of gases

        Tmin,Tmax = 150, 450 #temperature
        error = 1e9 #we set the error 100 0000
        Ngas = 6 #number of gases
        Fmin = np.linspace(0.7,1,30)



        for i in xrange(Ngas): #iterating through the gases
            print i, 'iteration'
            z = -1
            for absorbline in self.gases()[i]: #iterating through the absorptionline of gas i
                z +=1
                """
                this loop iterates through every gas and their absorption spectre,
                and makes a truth and false array, so this can be used
                to find the interesting lambdas which is lambdashift + lambda0.
                """

                I1 = -1; I2 = -1; I3 = -1
                discovered = np.logical_and(self.spectrum[:,0] > absorbline - 0.1,
                self.spectrum[:,0] < absorbline + 0.1)

                #shifted absorption lines measured from satelite
                ltow, lback = self.satlambda(absorbline)
                lambdas = np.linspace(lback,ltow,300)
                #The divitation in the measurments, sigma
                width1 = self.width(self.FWHM(absorbline,Tmax,self.gasmass()[i]))
                width2 = self.width(self.FWHM(absorbline,Tmin,self.gasmass()[i]))
                sigma = np.linspace(width1,width2, 30)

                #interesting wavelenght, flux, and noise from the files
                wl = self.spectrum[:,0][discovered]
                obsflux = self.spectrum[:,1][discovered]
                signoise = self.sigmanoise[:,1][discovered]

                #print np.shape(lambdas), np.shape(Fmin)

                for f in Fmin: #minimal flux
                    I1 += 1; I2 = -1
                    for s in sigma: #diviation measurred from real values
                        I2 += 1; I3 = -1
                        for wlength in lambdas: #wavelength observed
                            I3 += 1
                            fluxmodel = self.flux(wl,f,s,wlength)
                            FLUXsquared = np.sum((obsflux-fluxmodel)**2)
                            lesserror = FLUXsquared/np.sum(signoise**2)
                            if lesserror < error:
                                error = lesserror
                                index = [I1,I2,I3]
                                #print index

                #print index
                #print np.shape(self.flux(wl,Fmin[index[0]],sigma[index[1]],lambdas[index[2]]))
                self.plot_lines(wl,obsflux,
                self.flux(lambdas[index[2]],Fmin[index[0]],sigma[index[1]],wl),absorbline)



















a = atmosphericparticles('spectrumfile', 'sigmafile')
a.interesting_lines_observed()
