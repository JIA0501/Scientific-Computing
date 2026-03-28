from ngsolve import *
from netgen.occ import *
import numpy as np
import os

# -------------------------
# User parameters
# -------------------------
Re = 400          # change this: 100, 200, 400
t_end = 80.0      # physical time
dt = 0.01         # time step
order_u = 3
order_p = 2
maxh = 0.02       # mesh size (coarser = faster)

Uin = 1.0
D = 0.1           # cylinder diameter (R=0.05)
nu = Uin * D / Re

# probe location (same idea as FD: downstream centerline-ish)
probe_x = 0.25
probe_y = 0.20

outdir = "fem_output"
os.makedirs(outdir, exist_ok=True)
probe_path = os.path.join(outdir, f"probe_fem_Re{Re}.csv")

# -------------------------
# Geometry: channel + cylinder
# -------------------------
L = 2.2
H = 0.41
cx, cy, r = 0.2, 0.2, 0.05

channel = Rectangle(L, H).Face()
cyl = Circle((cx, cy), r).Face()
geo = channel - cyl

# Name boundaries (important!)
# We use a simple heuristic: we'll later apply BCs by coordinate regions.
mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=maxh))

# -------------------------
# Finite element spaces
# -------------------------
V = VectorH1(mesh, order=order_u, dirichlet=".*")  # we'll enforce BCs strongly via Set on boundaries
Q = H1(mesh, order=order_p)
X = FESpace([V, Q])

(u, p), (v, q) = X.TnT()

# -------------------------
# Inflow profile
# -------------------------
# parabolic inflow in x-direction
uin = CoefficientFunction((4*(y*(H-y))/(H*H), 0))

# -------------------------
# Boundary regions via indicator functions
# left boundary: x=0
# right boundary: x=L (outflow)
# walls: y=0 or y=H
# cylinder: distance to (cx,cy) ~ r
# -------------------------
eps = 1e-6
is_left  = IfPos(eps - x, 1, 0)
is_right = IfPos(x - (L - eps), 1, 0)
is_bottom = IfPos(eps - y, 1, 0)
is_top = IfPos(y - (H - eps), 1, 0)
is_cyl = IfPos(r + 1e-3 - sqrt((x-cx)**2 + (y-cy)**2), 1, 0)

# For strong BC setting we will use "definedon" with BoundaryRegion
# We'll build BoundaryRegion by selecting points through coefficient functions is not straightforward,
# so we take a pragmatic approach: Set on named boundaries is hard without labels in OCC.
# Instead: enforce all Dirichlet on velocity by projecting BC using a mask in variational form.
# This keeps the script robust for your report purposes.

# -------------------------
# Time stepping: semi-implicit (diffusion implicit, convection explicit)
# -------------------------
gfu = GridFunction(X)
vel = gfu.components[0]
pres = gfu.components[1]

# initial condition
vel.Set((0, 0))
pres.Set(0)

vel_prev = vel.vec.CreateVector()
vel_prev.data = vel.vec

# bilinear form (implicit part)
a = BilinearForm(X, symmetric=False)
a += nu * InnerProduct(grad(u), grad(v)) * dx
a += (-div(v) * p - div(u) * q) * dx
a += 1e-10 * p * q * dx  # tiny pressure stabilization to avoid singularity

# mass form for velocity
m = BilinearForm(X, symmetric=True)
m += InnerProduct(u, v) * dx
m.Assemble()

a.Assemble()

# linear form
f = LinearForm(X)

# helper: apply "Dirichlet-like" mask by penalization (robust, quick to run)
# impose u=uin at left, u=0 at walls+cylinder
pen = 1e8
a_pen = BilinearForm(X, symmetric=True)
a_pen += pen * (is_left + is_bottom + is_top + is_cyl) * InnerProduct(u, v) * ds
a_pen.Assemble()

# RHS contributions for inflow
f_in = LinearForm(X)
f_in += pen * is_left * InnerProduct(uin, v) * ds
f_in.Assemble()

# system matrix for implicit step
A = m.mat.CreateMatrix()
A.AsVector().data = m.mat.AsVector() + dt * (a.mat.AsVector() + a_pen.mat.AsVector())

inv = A.Inverse(X.FreeDofs(), inverse="sparsecholesky")

# probe output
with open(probe_path, "w") as fp:
    fp.write("t,u,v\n")

t = 0.0
Draw(Norm(vel), mesh, "vel_norm")

# time loop
while t < t_end - 0.5*dt:
    t += dt

    # explicit convection term using previous velocity
    # (u_prev · ∇) u_prev
    u_prev = GridFunction(V)
    u_prev.vec.data = vel_prev

    conv = LinearForm(X)
    conv += -InnerProduct(grad(u_prev)*u_prev, v) * dx  # move to RHS
    conv.Assemble()

    rhs = gfu.vec.CreateVector()
    rhs.data = m.mat * gfu.vec + dt * (conv.vec + f_in.vec)

    # solve
    gfu.vec.data = inv * rhs

    # update prev
    vel_prev.data = vel.vec

    # sample probe
    uv = vel(probe_x, probe_y)
    with open(probe_path, "a") as fp:
        fp.write(f"{t},{uv[0]},{uv[1]}\n")

    if int(t/dt) % 200 == 0:
        print(f"t={t:.2f} Re={Re}")

print("DONE:", probe_path)
