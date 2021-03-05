import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from sys import exit
from theory import rho as rhoTH
from theory import rho2 as rho2TH
from theory import rhoN as rhoNTH
from theory import rho2N as rho2NTH


x0 = 0
L = 10.
T = 1.
l = 1.
v0 = 5.

phi1 = 0.
phi2 = 135


phi1 = phi1*2*np.pi/360
phi2 = phi2*2*np.pi/360
print(phi1,phi2)

rho = np.mean(np.loadtxt("rho.dat"),0)*L
bins = np.loadtxt("rho_bins.dat")
#theta = np.mean(np.loadtxt("theta.dat"),0)
#theta_bins = np.loadtxt("theta_bins.dat")

x = np.linspace(x0,x0+L, 1000)
r = rhoTH(x,L,T,l,phi1,phi2,v0)
N = rhoNTH(x0,L,T,l,phi1,phi2,v0)
r /= N

r2 = rho2TH(x,L,T,l,phi1,phi2,v0)
N = rho2NTH(x0,L,T,l,phi1,phi2,v0)
r2 /= N

plt.scatter(bins,rho)
plt.plot(x,r)
#plt.plot(x,r2)

plt.show()
