import numpy as np
import matplotlib.pyplot as plt
from glob import glob

def plot_one(file, label):
    a = np.loadtxt(file)
    k = a[:,0]
    d = a[:,1]
    plt.plot(k, d, label=label)

plt.figure()

plot_one("delta_jacobi_N50.dat", "Jacobi")
plot_one("delta_gs_N50.dat", "Gauss-Seidel")

for f in sorted(glob("delta_sor_N50_w*.dat")):
    w = f.split("w")[-1].replace(".dat","")
    plot_one(f, f"SOR ω={w}")

plt.yscale("log")         # log-lin：y log, x linear
plt.xlabel("iteration k")
plt.ylabel("delta(k)")
plt.title("Convergence measure δ vs iteration (N=50)")
plt.legend()
plt.tight_layout()
plt.show()
