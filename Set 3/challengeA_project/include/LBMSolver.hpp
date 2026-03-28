#pragma once
#include "Config.hpp"
#include "Grid.hpp"
#include "Probe.hpp"
#include <vector>
#include <cstdint>

class LBMSolver {
public:
  explicit LBMSolver(const Config& cfg);
  bool run(); // true stable, false unstable

private:
  Config cfg_;
  Grid g_;

  std::vector<double> f_, f2_; // Q=9 per cell
  std::vector<double> rho_, ux_, uy_;
  std::vector<uint8_t> solid_;

  Probe probe_;

  void init();
  void macroscopic();
  void collide();
  void stream();
  void apply_bc();
  void write_out(int step);

  inline int base(int i,int j) const { return 9*(j*g_.nx+i); }
};
