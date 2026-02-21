#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>

int main() {
    //Initial value
    int N = 100;
    double dx = 1.0/N;
    double dt = 0.001;
    double c = 1.0;
    double r = std::pow(c*dt/dx,2);

    //Open output file
    std::ofstream out("wave_iii.txt");

    //Data vector based on time level:n+1,n,n-1
    std::vector<double> psi_prev(N+1,0.0);
    std::vector<double> psi_curr(N+1,0.0);
    std::vector<double> psi_next(N+1,0.0);


    //Initail psi at t=0
    // Ψ(x,0) = sin(2πx)
    // and boundary: Ψ(0,t)=Ψ(1,t)=0
    for (int i= 0; i<=N; i++){
        double x = i*dx;
        //psi_curr[i] = sin(2*M_PI*x); //do three sin waves: sin(2πx), sin(5πx), local sin(5πx)
        //psi_curr[i] = sin(5*M_PI*x); 
        if (x > 0.2 && x < 0.4) psi_curr[i] = std::sin(5 * M_PI * x);
        else psi_curr[i] = 0.0;
    }
    psi_curr[0] = 0.0;
    psi_curr[N] = 0.0;

    //Calculate psi at t = dt
    for (int i = 1; i < N; i++){
        psi_next[i] = psi_curr[i] + 0.5*r*(psi_curr[i + 1] - 2.0*psi_curr[i] + psi_curr[i - 1]);
    }
    psi_next[0] = 0.0;
    psi_next[N] = 0.0;


    //Shift time levels: pervious<- curremt, current <- next.
    psi_prev = psi_curr;
    psi_curr = psi_next;

    //TIme Stepping
    int steps = 5000;
    int save_steps = 500;

    for(int n = 1; n<=steps; n++){
        for(int i = 1; i < N; i++){
            psi_next[i] = 2.0*psi_curr[i] - psi_prev[i] + r*(psi_curr[i + 1] - 2.0*psi_curr[i] + psi_curr[i - 1]);
        }

        psi_next[0] = 0.0;
        psi_next[N] = 0.0;

        //Simple sanity check: print data per 200 steps
        if (n % 1000 == 0) {
            double max_psi = 0.0;
            for (double v : psi_next) max_psi = std::max(max_psi, std::abs(v));
            std::cout<< "Step:" << n
                        <<", t = " << (n+1) * dt 
                        <<", Psi[mid]=" << psi_next[N/2]
                        <<", maxMax Psi:" << max_psi
                        <<"\n";
        }

        //Save data per 500 steps
        if(n % save_steps == 0) {
            double t = (n+1) * dt;

            out << "# t = " << t << "\n";

            for (int i = 0; i <= N; i++) {
            double x = i * dx;
            out << x << " " << psi_next[i] << "\n";
        }

        out<< "\n";
    }
        psi_prev = psi_curr;
        psi_curr = psi_next;
    }

    out.close();
    return 0;
}

