import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from sys import exit



rho = np.mean(np.loadtxt("rho.dat"),0)
bins = np.loadtxt("rho_bins.dat")

fx = np.mean(np.loadtxt("fx.dat"),0)
fy = np.mean(np.loadtxt("fy.dat"),0)
fz = np.mean(np.loadtxt("fz.dat"),0)

px = np.mean(np.loadtxt("px.dat"),0)
py = np.mean(np.loadtxt("py.dat"),0)
pz = np.mean(np.loadtxt("pz.dat"),0)

Dr = 20.
f = 10.
dT = f*f/(6*Dr)
L = 10.
w = 2*np.pi/L
k = 1.
print(dT)
def r(x):
    return 1./np.sqrt(1+dT*np.sin(w*x)**2)

def px_th(x):
    I1 = 1/(1+dT)
    I2 = 1/(1+k**2)
    I3 = 1./(6*Dr)
    return -I1*I2*I3*f*w*np.cos(w*x)

def py_th(x):
    return k*px_th(x)

rr = r(bins)
norm = simps(rr,bins)
rr /= norm

norm = simps(rho,bins)
rho /= norm

plt.subplot(1,3,1)
plt.plot(bins,rho)
plt.plot(bins,rr)
plt.title("density")

plt.subplot(1,3,2)
plt.plot(bins,px)
plt.plot(bins,px_th(bins))
plt.title("px")


plt.subplot(1,3,3)
plt.plot(bins,py)
plt.plot(bins,py_th(bins))
plt.title("py")


plt.show()

exit()
plt.subplot(1,3,2)
plt.plot(bins,fx, label="fx")
plt.plot(bins,fy, label="fy")
plt.legend()
plt.title("flux")


