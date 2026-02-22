#include "iteration.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>


Laplace2D::Laplace2D(int N_)
    :N(N_), Nx(N_), Ny(N_ + 1){
    
    dx = 1.0/static_cast<double>(N);
    
    c.assign(Nx*Ny, 0.0);
    c_new.assign(Nx*Ny, 0.0);

    obj.assign(Nx*Ny, 0);
    enforce_bc(c);
    enforce_bc(c_new);

    
    
}

void Laplace2D::enforce_bc(std::vector<double>& a){
    for(int i=0; i<Nx; ++i) a[idx(i, 0)] = 0.0;
    for(int i=0; i<Nx; ++i) a[idx(i, Ny-1)] = 1.0;
}

void Laplace2D::add_rect(int i0, int i1, int j0, int j1)
{
    
    j0 = std::max(0, j0);
    j1 = std::min(Ny - 1, j1);

    for (int j = j0; j <= j1; ++j) {
        for (int i = i0; i <= i1; ++i) {
            int ii = (i % Nx + Nx) % Nx;  
            obj[idx(ii, j)] = 1;
        }
    }
}

void Laplace2D::enforce_sink(std::vector<double>& a)
{
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            if (is_obj(i, j)) a[idx(i, j)] = 0.0;
        }
    }
}


//Jacobi iterarion
int jacobi(Laplace2D& g, double eps, int max_iters, double& last_del, std::ostream* trace){
    int it = 0;

    for(; it<max_iters; ++it){
        double delta = 0.0;

        for (int j = 1; j <=g.Ny-2; ++j){
            for(int i=0; i<g.Nx; ++i){

                if (g.is_obj(i,j)) { g.c_new[g.idx(i,j)] = 0.0; continue; }

                int ip = (i+1)%g.Nx;
                int im = (i-1+g.Nx)%g.Nx;

                double newv = 0.25 * (
                    g.c[g.idx(ip,j)] + g.c[g.idx(im,j)] + g.c[g.idx(i,j+1)] + g.c[g.idx(i,j-1)]

                );

                double diff = std::abs(newv-g.c[g.idx(i,j)]);
                if(diff>delta) delta = diff;

                g.c_new[g.idx(i,j)] = newv;
            }
        }
        g.enforce_bc(g.c_new);
        g.enforce_sink(g.c_new);
        g.c.swap(g.c_new);

        last_del = delta;
        if (trace) (*trace) << it << " " << delta << "\n";
        if(delta < eps) break;
    }

    return it+1;
}

//Gauss-Seidel iteration
int gauss_seidel(Laplace2D& g, double eps, int max_iter, double& last_delta, std::ostream* trace)
{
    int it = 0;
    for(; it < max_iter; ++it){
        double delta = 0.0;

        for(int j=1; j<=g.Ny-2; ++j){
            for(int i=0; i<g.Nx; ++i){

                if (g.is_obj(i,j)) { g.c[g.idx(i,j)] = 0.0; continue; }

                int ip = (i+1) % g.Nx;
                int im = (i-1+g.Nx) % g.Nx;

                double old = g.c[g.idx(i,j)];
                double newv = 0.25 * (
                    g.c[g.idx(ip,j)] + g.c[g.idx(im,j)] +
                    g.c[g.idx(i,j+1)] + g.c[g.idx(i,j-1)]
                );

                double diff = std::abs(newv - old);
                if(diff > delta) delta = diff;

                g.c[g.idx(i,j)] = newv;
            }
        }

        g.enforce_bc(g.c);
        g.enforce_sink(g.c);
        last_delta = delta;
        if (trace) (*trace) << it << " " << delta << "\n";
        if(delta < eps) break;
    }
    return it+1;
}

//SOR iteration
int sor(Laplace2D& g, double omega, double eps, int max_iter, double& last_delta, std::ostream* trace)
{
    int it = 0;
    for(; it < max_iter; ++it){
        double delta = 0.0;

        for(int j=1; j<=g.Ny-2; ++j){
            for(int i=0; i<g.Nx; ++i){

                if (g.is_obj(i,j)) { g.c[g.idx(i,j)] = 0.0; continue; }

                int ip = (i+1) % g.Nx;
                int im = (i-1+g.Nx) % g.Nx;

                double old = g.c[g.idx(i,j)];
                double gs  = 0.25 * (
                    g.c[g.idx(ip,j)] + g.c[g.idx(im,j)] +
                    g.c[g.idx(i,j+1)] + g.c[g.idx(i,j-1)]
                );

                double newv = omega * gs + (1.0 - omega) * old;

                double diff = std::abs(newv - old);
                if(diff > delta) delta = diff;

                g.c[g.idx(i,j)] = newv;
            }
        }

        g.enforce_bc(g.c);
        g.enforce_sink(g.c);
        last_delta = delta;
        if (trace) (*trace) << it << " " << delta << "\n";
        if(delta < eps) break;
    }
    return it+1;
}

double max_err_vs_linear(const Laplace2D& g)
{
    int i_mid = g.Nx / 2;      
    double emax = 0.0;

    for(int j = 0; j < g.Ny; ++j){
        double y = j * g.dx;               
        double num = g.c[g.idx(i_mid, j)]; 
        double ana = y;                   
        emax = std::max(emax, std::abs(num - ana));
    }
    return emax;
}

void write_field_xyz(const Laplace2D& g, const std::string& filename)
{
    std::ofstream out(filename);
    out << std::setprecision(12);

    for(int j=0; j<g.Ny; ++j){
        double y = j * g.dx;
        for(int i=0; i<g.Nx; ++i){
            double x = i * g.dx;
            out << x << " " << y << " " << g.c[g.idx(i,j)] << "\n";
        }
    }
}
