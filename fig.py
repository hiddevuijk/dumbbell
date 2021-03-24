import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from sys import exit
from numerical import rho as rhoNUM
from numerical import rhoN as rhoNNUM

from theory import rho as rhoTH
from theory import rhoN as rhoNTH



x0 = 0
L = 10.
T = 1.
l = 1.
v0 = 20.

phi1 = 0.
phi2 = 1.57

sm = np.sin(phi1) - np.sin(phi2) 
sp = np.sin(phi1) + np.sin(phi2) 

cm = np.cos(phi1) - np.cos(phi2) 
cp = np.cos(phi1) + np.cos(phi2) 

#print(cm, cp, sm,sp)


rho = np.mean(np.loadtxt("rho.dat"),0)*L
bins = np.loadtxt("rho_bins.dat")
#theta = np.mean(np.loadtxt("theta.dat"),0)
#theta_bins = np.loadtxt("theta_bins.dat")

x = np.linspace(x0,x0+L, 1000)
r = rhoTH(x,L,T,l,phi1,phi2,v0)
N = rhoNTH(x0,L,T,l,phi1,phi2,v0)
r /= N/L


#rnum = np.asarray( [ rhoNUM(xi,L,T,l,phi1,phi2,v0) for xi in x] )
#N = rhoNNUM(L,T,l,phi1,phi2,v0)
#rnum /= N/L

plt.subplot(1,2,1)
plt.scatter(bins,rho*L)
plt.plot(x,r, color="red")
#plt.plot(x,rnum, color="black", linestyle=":")
#plt.plot(x,r2, color="blue")

plt.subplot(1,2,2)

#Jx = np.mean( np.loadtxt("fluxX.dat"), axis = 0)
#Jy = np.mean( np.loadtxt("fluxY.dat") , axis = 0)
#plt.scatter(bins, Jx*L, label="Jx")
#plt.scatter(bins, Jy*L, label="Jy")

px = np.mean( np.loadtxt("px.dat"), axis=0)
py = np.mean( np.loadtxt("py.dat"), axis=0)
plt.scatter(bins,px, label="px")
plt.scatter(bins,py, label="py")
plt.legend()

plt.tight_layout()
plt.show()
