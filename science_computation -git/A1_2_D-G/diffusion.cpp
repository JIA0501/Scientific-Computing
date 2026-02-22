#include "diffusion.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

Diffusion2D::Diffusion2D(int N_, double D_, double dt_)
: N(N_), Nx(N_), Ny(N_ + 1), D(D_), dt(dt_)
{
    dx = 1.0 / static_cast<double>(N);
    alpha = dt * D / (dx * dx);

    c.assign(Nx * Ny, 0.0);
    c_new.assign(Nx * Ny, 0.0);

    enforce_dirichlet_y(c);
    enforce_dirichlet_y(c_new);
}

int Diffusion2D::idx(int i, int j) const { return i + j * Nx; }

void Diffusion2D::enforce_dirichlet_y(std::vector<double>& a) {
    for (int i = 0; i < Nx; ++i) a[idx(i, 0)] = 0.0;
    for (int i = 0; i < Nx; ++i) a[idx(i, Ny - 1)] = 1.0;
}

void Diffusion2D::step() {
    for (int j = 1; j <= Ny - 2; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int ip = (i + 1) % Nx;
            int im = (i - 1 + Nx) % Nx;

            double cij = c[idx(i, j)];
            double lap = c[idx(ip, j)] + c[idx(im, j)]
                       + c[idx(i, j + 1)] + c[idx(i, j - 1)]
                       - 4.0 * cij;

            c_new[idx(i, j)] = cij + alpha * lap;
        }
    }
    enforce_dirichlet_y(c_new);
    c.swap(c_new);
}

double analytic(double y, double t, double D) {
    if (t <= 0.0) {
        std::cerr << "Error: t must be positive.\n";
        return 0.0;
    }
    double sum = 0.0;
    int M = 20;
    for (int i = 0; i < M; ++i) {
        double a = (1.0 - y + 2*i) / (2 * std::sqrt(D*t));
        double b = (1.0 + y + 2*i) / (2 * std::sqrt(D*t));
        sum += std::erfc(a) - std::erfc(b);
    }
    return sum;
}

void write_profile(const Diffusion2D& sim, double t_target, const std::string& filename) {
    std::ofstream out(filename);
    out << std::setprecision(12);

    int i_mid = sim.Nx / 2;
    for (int j = 0; j < sim.Ny; ++j) {
        double y = j * sim.dx;
        double num = sim.c[sim.idx(i_mid, j)];
        double ana = analytic(y, t_target, sim.D);
        out << y << " " << num << " " << ana << "\n";
    }
}

void write_field_xyz(const Diffusion2D& sim, const std::string& filename) {
    std::ofstream out(filename);
    out << std::setprecision(12);

    for (int j = 0; j < sim.Ny; ++j) {
        double y = j * sim.dx;
        for (int i = 0; i < sim.Nx; ++i) {
            double x = i * sim.dx;
            out << x << " " << y << " " << sim.c[sim.idx(i, j)] << "\n";
        }
    }
}
