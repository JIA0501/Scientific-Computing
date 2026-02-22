#include "iteration.hpp"
#include <iostream>
#include <fstream>

int main(){
    double eps = 1e-5;
    int max_iter = 5000000;

    std::ofstream out("optimal_omega.dat");
    out << "# N omega_best it_best\n";

    for(int N: {20, 30, 40, 50, 60, 80, 100}){
        double best_w = 1.0;
        int best_it = max_iter;
        double best_last = 0.0;

        for(double w = 1.0; w < 2.0; w += 0.02){
            Laplace2D g(N);
            double last=0.0;
            int it = sor(g, w, eps, max_iter, last, nullptr);
            if(it < best_it){
                best_it = it;
                best_w = w;
                best_last = last;
            }
        }

        std::cout << "N="<<N<<" best_w="<<best_w<<" best_it="<<best_it<<"\n";
        out << N << " " << best_w << " " << best_it << "\n";
    }
}
