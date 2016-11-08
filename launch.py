from __future__ import division
import numpy as np
from AST1100SolarSystem import AST1100SolarSystem
import scipy.interpolate as inter
import time
start = time.time()

system = AST1100SolarSystem(23078) #my seed
txtfile = ('part5.txt')
system.sendSatellite(txtfile)
