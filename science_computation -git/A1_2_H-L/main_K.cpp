#include "iteration.hpp"
#include <iostream>
#include <vector>
#include <string>

struct RunResult {
    std::string name;
    int it_best;
    double w_best;
};

RunResult run_one(const std::string& tag, Laplace2D& g, double eps, int max_iter)
{
    // Take a small-scale scan to observe the trend
    double best_w = 1.0;
    int best_it = max_iter;

    for(double w : {1.2,1.4,1.6,1.7,1.8,1.85,1.9}){
        Laplace2D tmp = g;            
        tmp.enforce_sink(tmp.c);

        double last=0;
        int it = sor(tmp, w, eps, max_iter, last, nullptr);
        if(it < best_it){
            best_it = it;
            best_w = w;
        }
    }

    // Use best_w run againly，Output the final field
    {
        double last=0;
        int it = sor(g, best_w, eps, max_iter, last, nullptr);
        std::cout << tag << ": best_w=" << best_w << " it=" << it << "\n";
        write_field_xyz(g, "K_" + tag + ".xyz");
    }

    return {tag, best_it, best_w};
}

int main(){
    int N = 50;
    double eps = 1e-5;
    int max_iter = 5000000;

    // A) baseline
    {
        Laplace2D g(N);
        run_one("baseline", g, eps, max_iter);
    }

    // B) one rectangle sink
    {
        Laplace2D g(N);
        g.add_rect(N/3, N/3+6, N/2, N/2+12);
        g.enforce_sink(g.c);
        run_one("rect1", g, eps, max_iter);
    }

    // C) two rectangles sink
    {
        Laplace2D g(N);
        g.add_rect(N/4, N/4+6, N/3, N/3+10);
        g.add_rect(3*N/5, 3*N/5+8, N/2, N/2+12);
        g.enforce_sink(g.c);
        run_one("rect2", g, eps, max_iter);
    }

    return 0;
}
