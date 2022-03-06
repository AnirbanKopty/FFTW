# TODO: To solve heat equation del/del t (theta) = kappa nabla^2 (theta)  (kappa = 1)
# Using Spectral Methods by FFTW

# We first go to the Fourier Space, thus nabla^2 = -k^2, then solve for t in that space and
# Finally, return back to real space to get the final solution

#? Here, we will be dealing with 2D array, but the fourier transformation will be on the space dimension only which is 1D

# Libraries
import numpy as np
import pyfftw
import matplotlib.pyplot as plt

# Dimensions
T = 10; L = 2*np.pi
Nt = 100; Nx = 8

dt = T/Nt; d


# Datatypes
theta = pyfftw.empty_aligned((), dtype='float64')