import numpy as np
import matplotlib.pyplot as plt
from sys import exit



rho = np.loadtxt("rho.dat")
bins = np.loadtxt("rho_bins.dat")

theta = np.loadtxt("theta.dat")

plt.subplot(1,2,1)
plt.imshow(rho,interpolation='none')
plt.colorbar()
plt.title("density")

plt.subplot(1,2,2)
plt.imshow(theta, interpolation='none')
plt.colorbar()
plt.title("theta")

plt.tight_layout()
plt.show()


