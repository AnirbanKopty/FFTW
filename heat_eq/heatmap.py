import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb

heat_sol = np.loadtxt("heat_sol_cmplx_ps.dat")

sb.heatmap(heat_sol)
plt.savefig("heat_sol.png")
# plt.show()