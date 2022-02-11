import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb

heat_sol = np.loadtxt("heat_sol.dat")

sb.heatmap(heat_sol)
# plt.show()
plt.savefig("heat_sol.png")