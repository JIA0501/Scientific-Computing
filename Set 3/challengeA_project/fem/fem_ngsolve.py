import argparse
import os
from ngsolve import *
from netgen.geom2d import SplineGeometry

def make_geometry(Lx, Ly, cx, cy, R):
    geo = SplineGeometry()
    p1 = geo.AppendPoint(0, 0)
    p2 = geo.AppendPoint(Lx, 0)
    p3 = geo.AppendPoint(Lx, Ly)
    p4 = geo.AppendPoint(0, Ly)
    geo.Append(["line", p1, p2], bc="bottom")
    geo.Append(["line", p2, p3], bc="outlet")
    geo.Append(["line", p3, p4], bc="top")
    geo.Append(["line", p4, p1], bc="inlet")
    geo.AddCircle((cx, cy), R, bc="cyl", leftdomain=0, rightdomain=1)
    return geo

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Re", type=float, default=100.0)
    ap.add_argument("--steps", type=int, default=2000)
    ap.add_argument("--dt", type=float, default=1e-3)
    ap.add_argument("--maxh", type=float, default=0.02)
    ap.add_argument("--outdir", type=str, default="output")
    args = ap.parse_args()

    Lx, Ly = 2.2, 0.41
    cx, cy, R = 0.2, 0.205, 0.05
    Uin = 1.0
    nu = Uin*(2*R)/args.Re

    geo = make_geometry(Lx, Ly, cx, cy, R)
    mesh = Mesh(geo.GenerateMesh(maxh=args.maxh))
    mesh.Curve(3)

    V = VectorH1(mesh, order=2, dirichlet="inlet|top|bottom|cyl")
    Q = H1(mesh, order=1)
    X = FESpace([V, Q])

    (u,p) = X.TrialFunction()
    (v,q) = X.TestFunction()

    y = ycoord(mesh)
    inflow = CoefficientFunction((4*Uin*y*(Ly-y)/(Ly*Ly), 0))

    gfu = GridFunction(X)
    gfu.components[0].Set(CoefficientFunction((0,0)))
    gfu.components[0].Set(inflow, definedon=mesh.Boundaries("inlet"))
    gfu.components[1].Set(0.0)

    dt = args.dt

    a = BilinearForm(X, symmetric=False)
    a += (nu*InnerProduct(grad(u), grad(v)) - div(v)*p - q*div(u)) * dx
    gamma = 1e-3
    a += gamma*div(u)*div(v) * dx
    a.Assemble()

    m = BilinearForm(X)
    m += InnerProduct(u, v) * dx
    m.Assemble()

    os.makedirs(args.outdir, exist_ok=True)
    vtk = VTKOutput(mesh, coefs=[gfu.components[0], gfu.components[1]],
                    names=["u","p"], filename=os.path.join(args.outdir,"fem"),
                    subdivision=2)

    free = X.FreeDofs()

    for it in range(args.steps):
        un = gfu.components[0]
        conv = InnerProduct(Grad(un)*un, v) * dx

        rhs = LinearForm(X)
        rhs += (1/dt)*InnerProduct(un, v)*dx
        rhs += (-conv)
        rhs.Assemble()

        A = (1/dt)*m.mat + a.mat
        inv = A.Inverse(free, inverse="umfpack")
        gfu.vec.data = inv * rhs.vec

        gfu.components[0].Set(CoefficientFunction((0,0)))
        gfu.components[0].Set(inflow, definedon=mesh.Boundaries("inlet"))

        if it % 50 == 0:
            print(f"[FEM] step {it} Re={args.Re} nu={nu}")
            vtk.Do()

if __name__ == "__main__":
    main()
