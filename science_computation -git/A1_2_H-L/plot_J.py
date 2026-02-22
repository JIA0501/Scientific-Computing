import numpy as np
import matplotlib.pyplot as plt

a = np.loadtxt("optimal_omega.dat")
N = a[:,0]
w = a[:,1]
it = a[:,2]

plt.figure()
plt.plot(N, w, marker="o")
plt.xlabel("N")
plt.ylabel("optimal ω")
plt.title("Optimal ω vs N (SOR)")
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(N, it, marker="o")
plt.xlabel("N")
plt.ylabel("iterations to converge")
plt.title("Iterations at optimal ω vs N")
plt.tight_layout()
plt.show()
