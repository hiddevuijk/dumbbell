import numpy as np
import matplotlib.pyplot as plt
from sys import exit



rho = np.loadtxt("rho.dat")
bins = np.loadtxt("rho_bins.dat")

theta = np.loadtxt("theta.dat")

Jx = np.loadtxt("fluxX.dat")
Jy = np.loadtxt("fluxY.dat")

plt.subplot(1,2,1)
plt.imshow(Jx,interpolation='none')
plt.colorbar()
plt.title("Jx")

plt.subplot(1,2,2)
plt.imshow(Jy, interpolation='none')
plt.colorbar()
plt.title("Jy")

plt.tight_layout()
plt.show()


