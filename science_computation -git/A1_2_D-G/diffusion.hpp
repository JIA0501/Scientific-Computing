#pragma once
#include <vector>
#include <string>

struct Diffusion2D {
    int N;
    int Nx, Ny;
    double D;
    double dx, dt;
    double alpha;
    std::vector<double> c, c_new;

    Diffusion2D(int N_, double D_, double dt_);
    int idx(int i, int j) const;
    void enforce_dirichlet_y(std::vector<double>& a);
    void step();
};

double analytic(double y, double t, double D);

void write_profile(const Diffusion2D& sim, double t_target, const std::string& filename);


void write_field_xyz(const Diffusion2D& sim, const std::string& filename);
