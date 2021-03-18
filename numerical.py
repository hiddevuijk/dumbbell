import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from scipy.integrate import quad



def I(x, L, T, l, phi1, phi2, v0):
    f = v0*(1 + np.sin( 2*np.pi*x/L + 3*np.pi/2 ) )
    df = v0* np.cos( 2*np.pi*x/L + 3*np.pi/2 ) *2*np.pi/L

    sm = np.sin(phi1) - np.sin(phi2) 
    sp = np.sin(phi1) + np.sin(phi2) 

    cm = np.cos(phi1) - np.cos(phi2) 
    cp = np.cos(phi1) + np.cos(phi2) 

    a = l*f/(2*T)
    da = l*df/(2*T)

    num = cm - (cp*cp+2*sp*sp + sm*sp*cp*a)*a/(1+a*a*sm*sm)
    num *= l*da
    den = 2+a*a*(cp*cp+sp*sp)/(1+a*a*sm*sm)

    return num/den


def rho(x, L, T, l, phi1, phi2, v0):
    return np.exp( quad( lambda y: I(y,L,T,l,phi1,phi2,v0), 0,x)[0] )

def rhoN(L, T, l, phi1, phi2, v0):
    return quad( lambda y: rho(y, L,T,l,phi1, phi2,v0), 0,L)[0]



