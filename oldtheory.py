import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt


    

def rho(x, L, T, l, phi1, phi2, v0):

    f = v0*(1 + np.sin( 2*np.pi*x/L + 3*np.pi/2 ) )

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
    

def rho2(x, L, T, l, phi1, phi2, v0):

    f = v0*(1 + np.sin( 2*np.pi*x/L + 3*np.pi/2) )

    d = np.cos(phi1) - np.cos(phi2)

    A = l*d*f/(4*T)
    r = np.exp(A)
    return r


def rhoN(x0, L, T, l, phi1, phi2, v0):
    return quad( lambda x: rho(x,L,T,l,phi1,phi2,v0), x0, x0+L)[0]
 
def rho2N(x0, L, T, l, phi1, phi2, v0):
    return quad( lambda x: rho2(x,L,T,l,phi1,phi2,v0), x0, x0+L)[0]
       
