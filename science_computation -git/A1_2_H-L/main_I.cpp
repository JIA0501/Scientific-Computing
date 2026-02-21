#include "iteration.hpp"
#include <fstream>
#include <iostream>

int main(){
    int N = 50;
    double eps = 1e-5;
    int max_iter = 2000000;

    // Jacobi
    {
        Laplace2D g(N);
        std::ofstream tr("delta_jacobi_N50.dat");
        double last=0.0;
        int it = jacobi(g, eps, max_iter, last, &tr);
        std::cout << "Jacobi it=" << it << " last_delta=" << last << "\n";
    }

    // GS
    {
        Laplace2D g(N);
        std::ofstream tr("delta_gs_N50.dat");
        double last=0.0;
        int it = gauss_seidel(g, eps, max_iter, last, &tr);
        std::cout << "GS it=" << it << " last_delta=" << last << "\n";
    }

    // SOR
    for(double omega: {1.2, 1.5, 1.7, 1.8, 1.9}){
        Laplace2D g(N);
        std::ofstream tr("delta_sor_N50_w" + std::to_string(omega) + ".dat");
        double last=0.0;
        int it = sor(g, omega, eps, max_iter, last, &tr);
        std::cout << "SOR w="<<omega<<" it="<<it<<" last_delta="<<last<<"\n";
    }
}
