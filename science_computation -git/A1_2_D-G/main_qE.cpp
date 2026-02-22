#include "diffusion.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

int main(){
    int N = 100;
    double D = 1.0;

    double dx = 1.0 / N;
    double dt_max = dx*dx / (4.0*D);
    double dt = 0.2 * dt_max;

    Diffusion2D sim(N, D, dt);

    std::vector<double> times = {0.001, 0.01, 0.1, 1.0};
    double current_time = 0.0;

    for(double t_target : times) {

        while(current_time < t_target){
            sim.step();
            current_time += sim.dt;
        }

        std::ostringstream oss;
        oss << "profile_t_" << std::fixed << std::setprecision(6) << t_target << ".dat";
        write_profile(sim, t_target, oss.str());

        std::cout << "Wrote " << oss.str() << "\n";
    }

    return 0;
}
