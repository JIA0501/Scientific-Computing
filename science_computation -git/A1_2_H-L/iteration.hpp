#pragma once
#include <vector>
#include <iosfwd>



struct Laplace2D {
        int N, Nx, Ny;
        double dx;
        std::vector<double> c, c_new;

        Laplace2D(int N_);
        inline int idx(int i, int j) const {return i+j*Nx;};
        void enforce_bc(std::vector<double>& a);
        void add_rect(int i0,int i1,int j0,int j1);
        void enforce_sink(std::vector<double>& a);

        std::vector<int> obj; // 0=No stuff, 1=stuffsink)
        bool is_obj(int i,int j) const { return obj[idx(i,j)] != 0; }

};

int jacobi(Laplace2D& g, double eps, int max_iters, double& last_del, std::ostream* trace=nullptr);
int gauss_seidel(Laplace2D& g, double eps, int max_iter, double& last_delta, std::ostream* trace=nullptr);
int sor(Laplace2D& g, double omega, double eps, int max_iter, double& last_delta, std::ostream* trace=nullptr);

double max_err_vs_linear(const Laplace2D& g);

// K: write field for plotting
void write_field_xyz(const Laplace2D& g, const std::string& filename);