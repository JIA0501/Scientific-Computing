import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

RUN_STEP = 5

# STEP 0: PARAMETERS

# Domain size
Lx, Ly = 10.0, 8.0

# Grid spacing
h = 0.1   # start coarse; later maybe 0.05

Nx = int(Lx / h) + 1
Ny = int(Ly / h) + 1

x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)
X, Y = np.meshgrid(x, y, indexing="xy")

# Physical constants
c0 = 3e8
f_eff = 0.8e9
k0 = 2 * np.pi * f_eff / c0

# Materials
n_air = 1.0 + 0.0j
n_wall = 2.5 + 0.5j

# Source parameters
A_src = 1e4
sigma = 0.2

# Measurement points
measurement_points = np.array([
    [1.0, 5.0],   # Living room
    [2.0, 1.0],   # Kitchen
    [9.0, 1.0],   # Bathroom
    [9.0, 7.0],   # Bedroom 1
])

measurement_radius = 0.05
router_exclusion = 0.5

# STEP 1: GEOMETRY / WALL MASK
def build_wall_mask():
    wall = np.zeros((Ny, Nx), dtype=bool)

    def rect(x0, x1, y0, y1):
        ix0 = max(0, int(round(x0 / h)))
        ix1 = min(Nx - 1, int(round(x1 / h)))
        iy0 = max(0, int(round(y0 / h)))
        iy1 = min(Ny - 1, int(round(y1 / h)))
        wall[iy0:iy1+1, ix0:ix1+1] = True

    t = 0.15

    # Outer walls
    rect(0.0, Lx, 0.0, t)
    rect(0.0, Lx, Ly - t, Ly)
    rect(0.0, t, 0.0, Ly)
    rect(Lx - t, Lx, 0.0, Ly)

    # Internal walls approximating the figure
    rect(0.0, 3.0, 3.0 - t/2, 3.0 + t/2)
    rect(4.0, 6.0, 3.0 - t/2, 3.0 + t/2)
    rect(7.0, 10.0, 3.0 - t/2, 3.0 + t/2)

    rect(6.0 - t/2, 6.0 + t/2, 3.0, 8.0)

    rect(7.0 - t/2, 7.0 + t/2, 2.5, 3.0)
    rect(2.5 - t/2, 2.5 + t/2, 0.0, 2.0)
    rect(7.0 - t/2, 7.0 + t/2, 0.0, 1.5)

    return wall


wall_mask = build_wall_mask()
air_mask = ~wall_mask
k_field = np.where(wall_mask, n_wall * k0, n_air * k0).astype(np.complex128)

def plot_wall_mask():
    plt.figure(figsize=(8, 6))
    plt.imshow(wall_mask, origin="lower", extent=[0, Lx, 0, Ly], aspect="equal")
    plt.colorbar(label="wall")
    plt.scatter(
        measurement_points[:, 0], measurement_points[:, 1],
        marker="o", s=70, facecolors="none", edgecolors="red", linewidths=2,
        label="measurement points"
    )
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.title("Wall mask")
    plt.legend()
    plt.tight_layout()
    plt.show()

# STEP 2: GAUSSIAN SOURCE
def gaussian_source(xr, yr):
    r2 = (X - xr) ** 2 + (Y - yr) ** 2
    return A_src * np.exp(-r2 / (2 * sigma ** 2))


def plot_source(xr, yr):
    src = gaussian_source(xr, yr)
    plt.figure(figsize=(8, 6))
    plt.imshow(src, origin="lower", extent=[0, Lx, 0, Ly], aspect="equal")
    plt.colorbar(label="source strength")
    plt.scatter([xr], [yr], marker="*", s=180, c="red", label="router")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.title(f"Step 2: Gaussian source at ({xr:.1f}, {yr:.1f})")
    plt.legend()
    plt.tight_layout()
    plt.show()

# STEP 3: HELMHOLTZ SOLVER
def idx(i, j):
    return j * Nx + i


def solve_helmholtz(xr, yr):
    """
    First working version:
        Δu + k(x,y)^2 u = f
    with simple Dirichlet-like boundary treatment on outer boundary.
    Later we can replace this with a better Robin/absorbing BC.
    """
    Ntot = Nx * Ny
    A = lil_matrix((Ntot, Ntot), dtype=np.complex128)
    b = np.zeros(Ntot, dtype=np.complex128)

    fsrc = gaussian_source(xr, yr)

    for j in range(Ny):
        for i in range(Nx):
            p = idx(i, j)
            kij = k_field[j, i]

            # Boundary nodes: simple placeholder treatment
            if i == 0 or i == Nx - 1 or j == 0 or j == Ny - 1:
                A[p, p] = 1.0
                b[p] = 0.0
                continue

            # Interior 5-point stencil
            A[p, idx(i, j)]     = -4.0 / h**2 + kij**2
            A[p, idx(i + 1, j)] =  1.0 / h**2
            A[p, idx(i - 1, j)] =  1.0 / h**2
            A[p, idx(i, j + 1)] =  1.0 / h**2
            A[p, idx(i, j - 1)] =  1.0 / h**2

            b[p] = fsrc[j, i]

    A = csr_matrix(A)
    u = spsolve(A, b)
    return u.reshape((Ny, Nx))


def plot_solution(u, xr, yr, title="WiFi field"):
    plt.figure(figsize=(9, 6))
    signal_db = 20 * np.log10(np.abs(u) + 1e-12)

    plt.imshow(signal_db, origin="lower", extent=[0, Lx, 0, Ly], aspect="equal")
    plt.colorbar(label="20 log10 |u|")

    wy, wx = np.where(wall_mask)
    plt.scatter(wx * h, wy * h, s=3, c="black", label="walls")

    plt.scatter(
        measurement_points[:, 0], measurement_points[:, 1],
        marker="o", s=70, facecolors="none", edgecolors="white", linewidths=2,
        label="measurement points"
    )

    plt.scatter([xr], [yr], marker="*", s=180, c="red", label="router")

    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.title(title)
    plt.legend(loc="upper right")
    plt.tight_layout()
    plt.show()

# STEP 4: SCORING
def measure_signal(u, xm, ym):
    r = 0.05  # 5 cm

    mask = (X - xm)**2 + (Y - ym)**2 <= r**2
    values = np.abs(u[mask])

    if len(values) == 0:
        return 0.0

    return np.mean(values)


def score_field(u):
    vals = []
    for px, py in measurement_points:
        vals.append(measure_signal(u, px, py))
    return np.sum(vals), vals

def distance_to_nearest_wall(xr, yr):
    wy, wx = np.where(wall_mask)
    wall_x = wx * h
    wall_y = wy * h
    d = np.sqrt((wall_x - xr)**2 + (wall_y - yr)**2)
    return np.min(d)

def valid_router_position(xr, yr):
    wall_clearance = 0.6

    # Keep away from outer boundary
    if xr < wall_clearance or xr > Lx - wall_clearance:
        return False
    if yr < wall_clearance or yr > Ly - wall_clearance:
        return False

    i = int(round(xr / h))
    j = int(round(yr / h))

    i = min(max(i, 0), Nx - 1)
    j = min(max(j, 0), Ny - 1)

    # Cannot be inside a wall
    if wall_mask[j, i]:
        return False

    # Keep away from measurement points
    d = np.sqrt((measurement_points[:, 0] - xr) ** 2 +
                (measurement_points[:, 1] - yr) ** 2)
    if np.any(d < router_exclusion):
        return False

    # NEW: keep away from walls
    if distance_to_nearest_wall(xr, yr) <= wall_clearance:
        return False

    return True

def score_router_position(xr, yr):
    u = solve_helmholtz(xr, yr)
    total, vals = score_field(u)
    return total, vals, u

#Step 5
def coarse_search(step=1.0):
    xs = np.arange(0.5, Lx - 0.49, step)
    ys = np.arange(0.5, Ly - 0.49, step)

    results = []

    for xr in xs:
        for yr in ys:
            if not valid_router_position(xr, yr):
                continue

            print(f"Testing ({xr:.2f}, {yr:.2f})...")

            u = solve_helmholtz(xr, yr)
            score, vals = score_field(u)

            results.append((score, xr, yr, vals, u))

    results.sort(key=lambda t: t[0], reverse=True)
    return results

def refine_search(x_center, y_center, half_width=1.0, step=0.5):
    xs = np.arange(x_center - half_width, x_center + half_width + 1e-12, step)
    ys = np.arange(y_center - half_width, y_center + half_width + 1e-12, step)

    results = []

    for xr in xs:
        for yr in ys:
            if not valid_router_position(xr, yr):
                continue

            print(f"Refining ({xr:.2f}, {yr:.2f})...")
            u = solve_helmholtz(xr, yr)
            score, vals = score_field(u)
            results.append((score, xr, yr, vals, u))

    results.sort(key=lambda t: t[0], reverse=True)
    return results

def plot_score_map(results):
    xs = [r[1] for r in results]
    ys = [r[2] for r in results]
    ss = [r[0] for r in results]

    plt.figure(figsize=(8, 6))
    plt.scatter(xs, ys, c=ss, s=80)
    plt.colorbar(label="score")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.title("Router candidate scores")
    plt.tight_layout()
    plt.show()
# MAIN
if __name__ == "__main__":
    print("RUN_STEP =", RUN_STEP)

    xr_test, yr_test = 5.0, 4.0

    if RUN_STEP >= 1:
        plot_wall_mask()

    if RUN_STEP >= 2:
        plot_source(xr_test, yr_test)

    print("Solving Helmholtz system for one test router position...")
    u_test = solve_helmholtz(xr_test, yr_test)
    print("Solved. Field shape =", u_test.shape)

    if RUN_STEP >= 3:
        plot_solution(u_test, xr_test, yr_test, "Step 3: Router at (5.0, 4.0)")

    if valid_router_position(xr_test, yr_test):
        total_score, point_scores = score_field(u_test)
        print("Finished Step 4: scoring.")
        print("Total score =", total_score)
        print("Per measurement point =", point_scores)
    else:
        print("Test router position is invalid.")

    print("\nRunning coarse search...")
    results = coarse_search(step=1.5)

    best_score, best_x, best_y, best_vals, best_u = results[0]

    print("\nBEST ROUTER POSITION")
    print("x =", best_x)
    print("y =", best_y)
    print("score =", best_score)
    print("per point =", best_vals)

    plot_solution(best_u, best_x, best_y, "Best router position")
    plot_score_map(results)

    print("\nRunning local refinement...")
    refined = refine_search(best_x, best_y, half_width=1.0, step=0.5)

    best_score, best_x, best_y, best_vals, best_u = refined[0]

    print("\nREFINED BEST ROUTER POSITION")
    print("x =", best_x)
    print("y =", best_y)
    print("score =", best_score)
    print("per point =", best_vals)

    plot_solution(best_u, best_x, best_y, "Refined best router position")
    plot_score_map(refined)