from matplotlib import projections
import numpy as np
import matplotlib.pyplot as plt

# To plot 3D
ax = plt.axes(projection='3d')      # Creates the 3D axes system

# tpts = np.arange(0,1000,1)
# xpts = np.linspace(0,2*np.pi,32)
x = np.linspace(-6,6,30)
y = np.linspace(-6,6,30)
# Creating x-y plane where the matrix values will be the z
X, Y = np.meshgrid(x,y)

def f(x,y):
    return np.sin(np.sqrt(x**2 + y**2))

Z = f(X,Y)

# ax.plot3D(X,Y,Z)
# plot3D is when, X,Y,Z are 1D arrays,
# while contour is when X,Y,Z are 2D arrays
ax.contour(X,Y,Z, 50)

print(np.shape(Z))

plt.show()