import numpy as np
import matplotlib.pyplot as plt

f_clean = np.loadtxt("f_clean.dat")
f_noisy = np.loadtxt("f_noisy.dat")
f_hat = np.loadtxt("f_hat.dat")
f_clean_hat = np.loadtxt("f_clean_hat.dat")
f_denoised = np.loadtxt("f_denoised.dat")

plt.rcParams['figure.figsize'] = [20, 5] # For all the plots

plt.plot(f_clean[:,0], f_clean[:,1], label='f_clean')
plt.plot(f_noisy[:,0], f_noisy[:,1], label='f_noisy')
plt.legend()

plt.figure()
plt.plot(f_hat[:,0], f_hat[:,1], label='f_hat')
plt.plot(f_clean_hat[:,0], f_clean_hat[:,1], label='f_clean_hat')
# plt.ylim(0,10)
plt.legend()

plt.figure()
plt.plot(f_denoised[:,0], f_denoised[:,1])
plt.ylim(-5,5)

plt.show()