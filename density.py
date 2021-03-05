import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt


    

def rho(x, L, T, l, phi1, phi2, v0):

    f = v0*(1 + np.sin( 2*np.pi*x/L) )

    c = 1 + np.cos(phi1-phi2)
    d = np.cos(phi1) - np.cos(phi2)
    s = np.sin(phi1)+ np.sin(phi2)

    A = 4*T*T + l*l*c*(f**2)

    B = -1/2 - s*s/(4*c)

    C = d/(2*np.sqrt(c) )
    C *= np.arctan( l*np.sqrt(c)*f/(2*T) )

    r = A**B
    r*= np.exp(C)
    return r

def rhoN(x0, L, T, l, phi1, phi2, v0):
    return quad( lambda x: rho(x,L,T,l,phi1,phi2,v0), x0, x0+L)[0]

x0 = 0
L = 10.
T = 1.
l = 1.
phi1 = 0
phi2 = 0
v0 = 1.

N = rhoN(x0,L,T,l,phi1,phi2,v0)
x = np.linspace(x0,L,1000)
y = rho(x,L,T,l,phi1,phi2,v0)/N

plt.plot(x,y)
plt.show()
    
