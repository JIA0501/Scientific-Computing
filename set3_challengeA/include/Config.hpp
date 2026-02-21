#pragma once
#include <string>

struct OutputCfg {
  int every = 200;
  std::string dir = "output";
};

struct ProbeCfg {
  int every = 10;
  double x = -1.0;
  double y = -1.0;
  std::string prefix = "probe";
};

struct SweepCfg {
  bool enabled = false;
  double startRe = 50.0;
  double growth = 1.25;
  double tol = 10.0;
  int maxRuns = 50;
  int warmup_steps = 5000;
};

struct Config {
  std::string method = "fd"; // fd | lbm | fem

  // geometry (typical cylinder-in-channel benchmark)
  double Lx = 2.2;
  double Ly = 0.41;
  double cx = 0.2;
  double cy = 0.205;
  double R  = 0.05;

  // discretization
  int nx = 800;
  int ny = 200;

  // physics
  double Re = 100.0;
  double Uin = 1.0;
  double rho = 1.0;

  // computed/solver-specific
  double nu  = 0.0;   // FD viscosity (computed)
  double dt  = 0.0;   // FD timestep (auto if 0)

  // time
  int steps = 20000;

  // FD params
  int poissonIters = 120;
  double poissonTol = 1e-6;

  // LBM params
  double tau = 0.0;    // auto if 0
  double U_lbm = 0.08; // inlet velocity in lattice units (keep small)

  // stability thresholds
  double speed_cap_factor = 20.0; // break if max speed > factor * Uin (FD) or factor * U_lbm (LBM)

  // output
  OutputCfg out;
  ProbeCfg probe;
  SweepCfg sweep;

  // threads
  int threads = 0; // 0 => OpenMP default
};
