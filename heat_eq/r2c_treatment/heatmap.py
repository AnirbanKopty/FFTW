from matplotlib import projections
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb

heat_sol = np.loadtxt("heat_sol.dat")
# heat_sol[t_pts,x_pts]

sb.heatmap(heat_sol)
plt.savefig("heat_sol.png")

# # To plot 3D (contour plot)
# plt.figure()
# ax = plt.axes(projection='3d')      # Creates the 3D axes system

# # tpts = np.arange(0,1000,1)
# # xpts = np.linspace(0,2*np.pi,32)
t = list(range(heat_sol.shape[0]))
x = list(range(heat_sol.shape[1]))
# # Creating x-y plane where the matrix values will be the z
# X, T = np.meshgrid(x,t)

# ax.contour(X,T,heat_sol,50)
# plt.savefig("heat_sol_3D.png")

# to plot waterfall plot (check the code)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
# plt.set_cmap('jet_r')
heat_plot = heat_sol[0:-1:100,:]    # choosing some time sections for plotting
for j in range(heat_plot.shape[0]):
    ys = j*np.ones(len(x))
    ax.plot(x,ys,heat_plot[j,:])#,color=cm.jet(j*20))

plt.savefig("heat_sol_waterfall.png")

plt.show()