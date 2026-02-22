#include "iteration.hpp"
#include <iostream>

int main() {
    int N = 50;
    double eps = 1e-5;
    int max_iter = 2000000;

    // Jacobi
    {
        Laplace2D g(N);
        double delta = 0.0;
        int it = jacobi(g, eps, max_iter, delta);
        std::cout << "Jacobi: it=" << it
                  << " delta=" << delta
                  << " max_err=" << max_err_vs_linear(g) << "\n";
    }

    // Gauss-Seidel
    {
        Laplace2D g(N);
        double delta = 0.0;
        int it = gauss_seidel(g, eps, max_iter, delta);
        std::cout << "Gauss-Seidel: it=" << it
                  << " delta=" << delta
                  << " max_err=" << max_err_vs_linear(g) << "\n";
    }

    // SOR
    {
        Laplace2D g(N);
        double delta = 0.0;
        double omega = 1.8;
        int it = sor(g, omega, eps, max_iter, delta);
        std::cout << "SOR(omega=" << omega << "): it=" << it
                  << " delta=" << delta
                  << " max_err=" << max_err_vs_linear(g) << "\n";
    }

    return 0;
}
