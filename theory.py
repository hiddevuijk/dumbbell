import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt


    

def rho(x, L, T, l, phi1, phi2, v0):

    f = v0*(1 + np.sin( 2*np.pi*x/L + 3*np.pi/2 ) )
    f *= l/(2*T)
    cp = np.cos(phi1) + np.cos(phi2)
    cm = np.cos(phi1) - np.cos(phi2)
    sp = np.sin(phi1) + np.sin(phi2)
    sm = np.sin(phi1) - np.sin(phi2)

    a = cm
    b = cp*cp+2*sp*sp
    c = cm*sm*sm-sm*sp*cp
    d = (2*sm*sm+cp*cp+sp*sp)/2

    A = (1 + d*f*f)**( -b/(4*d) )
    B = c*f/(2*d)
    C = (a*d-c)/(2*(d**(3./2)))
    C *= np.arctan( np.sqrt(d)*f)

    return A*np.exp( B+C)

def rhoN(x0, L, T, l, phi1, phi2, v0):
    return quad( lambda x: rho(x,L,T,l,phi1,phi2,v0), x0, x0+L)[0]
 
