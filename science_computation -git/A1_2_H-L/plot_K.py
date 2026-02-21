import numpy as np
import matplotlib.pyplot as plt
from glob import glob

N = 50
Nx, Ny = N, N+1

for fn in sorted(glob("K_*.xyz")):
    data = np.loadtxt(fn)
    C = data[:,2].reshape(Ny, Nx)

    plt.figure()
    im = plt.imshow(C, origin="lower", extent=[0,1,0,1], aspect="auto")
    plt.colorbar(im, label="c")
    plt.title(fn)
    plt.xlabel("x"); plt.ylabel("y")
    plt.tight_layout()
    plt.show()

