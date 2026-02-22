// wave_leapfrog.cpp
// Compile (C++20):
//   g++ -O2 -std=c++20 wave_leapfrog.cpp -o wave_leapfrog
// Run:
//   ./wave_leapfrog
//
// Output: wave_leapfrog.txt
//
// Leapfrog (staggered velocity):
//   v^{n+1/2} = v^{n-1/2} + dt * c^2 * Dxx(psi^n)
//   psi^{n+1} = psi^n + dt * v^{n+1/2}
// Bootstrap with v(x,0)=psi_t(x,0)=0:
//   v^{1/2} = v(0) + (dt/2)*c^2*Dxx(psi^0)

#include <cmath>
#include <numbers>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>

// Choose initial condition:
// 1: psi(x,0)=sin(2*pi*x)
// 2: psi(x,0)=sin(5*pi*x)
// 3: psi(x,0)=sin(5*pi*x) for 0.2<x<0.4 else 0
static int IC_TYPE = 3;

int main() {
    // Parameters
    const int N = 100;          // number of intervals => N+1 grid points
    const double L = 1.0;
    const double dx = L / N;
    const double dt = 0.001;
    const double c  = 1.0;

    const int steps = 5000;
    const int save_steps = 500;

    // Output file
    std::ofstream out("wave_leapfrog.txt");
    if (!out) {
        std::cerr << "ERROR: cannot open output file.\n";
        return 1;
    }

    // State: displacement psi^n at grid points i=0..N
    std::vector<double> psi(N + 1, 0.0);

    // Velocity at half-step: v^{n+1/2}
    std::vector<double> v_half(N + 1, 0.0);

    // ---- helper: discrete Laplacian in 1D ----
    auto laplacian = [&](const std::vector<double>& u, int i) -> double {
        return (u[i + 1] - 2.0 * u[i] + u[i - 1]) / (dx * dx);
    };

    // -------- Initial condition for psi(x,0) --------
    for (int i = 0; i <= N; i++) {
        const double x = i * dx;

        if (IC_TYPE == 1) {
            psi[i] = std::sin(2.0 * std::numbers::pi * x);
        } else if (IC_TYPE == 2) {
            psi[i] = std::sin(5.0 * std::numbers::pi * x);
        } else { // IC_TYPE == 3
            if (x > 0.2 && x < 0.4) psi[i] = std::sin(5.0 * std::numbers::pi * x);
            else psi[i] = 0.0;
        }
    }
    psi[0] = 0.0;
    psi[N] = 0.0;

    // -------- Bootstrap v^{1/2} from v(0)=0 --------
    for (int i = 1; i < N; i++) {
        v_half[i] = 0.5 * dt * (c * c) * laplacian(psi, i);
    }
    v_half[0] = 0.0;
    v_half[N] = 0.0;

    // Save initial state at t=0
    out << "# t = 0\n";
    for (int i = 0; i <= N; i++) {
        out << (i * dx) << " " << psi[i] << "\n";
    }
    out << "\n";

    // -------- Time stepping --------
    // Loop advances psi first (using v_half = v^{n+1/2}), then advances v_half to the next half-step.
    for (int n = 0; n < steps; n++) {
        const double t_next = (n + 1) * dt;

        // 1) Update displacement: psi^{n+1} = psi^n + dt * v^{n+1/2}
        for (int i = 1; i < N; i++) {
            psi[i] += dt * v_half[i];
        }
        psi[0] = 0.0;
        psi[N] = 0.0;

        // Diagnostics (every 1000 steps)
        if ((n + 1) % 1000 == 0) {
            double max_psi = 0.0;
            for (double val : psi) max_psi = std::max(max_psi, std::abs(val));
            std::cout << "Step: " << (n + 1)
                      << ", t = " << t_next
                      << ", Psi[mid] = " << psi[N / 2]
                      << ", max|Psi| = " << max_psi
                      << "\n";
        }

        // Save profile (every 500 steps)
        if ((n + 1) % save_steps == 0) {
            out << "# t = " << t_next << "\n";
            for (int i = 0; i <= N; i++) {
                out << (i * dx) << " " << psi[i] << "\n";
            }
            out << "\n";
        }

        // 2) Update velocity half-step:
        // v^{n+3/2} = v^{n+1/2} + dt * c^2 * Dxx psi^{n+1}
        for (int i = 1; i < N; i++) {
            v_half[i] += dt * (c * c) * laplacian(psi, i);
        }
        v_half[0] = 0.0;
        v_half[N] = 0.0;
    }

    out.close();
    return 0;
}