import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from sys import exit

from theory import rho as rhoTH
from theory import rhoN as rhoNTH


Nrho = 100
x0 = 0
L = 100.
T = 1.
l = 1.
v0 = 2.

phi1 = 0.
phi2 = 1.57

sm = np.sin(phi1) - np.sin(phi2) 
sp = np.sin(phi1) + np.sin(phi2) 

cm = np.cos(phi1) - np.cos(phi2) 
cp = np.cos(phi1) + np.cos(phi2) 

bins = np.loadtxt("rho_bins.dat")
for rhoi in range(0,Nrho,10):

    rho = np.loadtxt("results/rho{}.dat".format(rhoi))
    plt.scatter(bins,rho*L, label=rhoi)


# theory

x = np.linspace(x0,x0+L, 1000)
r = rhoTH(x,L,T,l,phi1,phi2,v0)
N = rhoNTH(x0,L,T,l,phi1,phi2,v0)
r /= N/L


plt.plot(x,r, color="red")

plt.legend()
plt.tight_layout()
plt.show()
