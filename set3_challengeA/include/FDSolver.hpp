#pragma once
#include "Config.hpp"
#include "Grid.hpp"
#include "Field.hpp"
#include "Probe.hpp"
#include <vector>
#include <cstdint>

class FDSolver {
public:
  explicit FDSolver(const Config& cfg);
  bool run(); // true stable, false unstable

private:
  Config cfg_;
  Grid g_;
  VectorField2 vel_, vel_tmp_;
  ScalarField p_, div_;
  std::vector<uint8_t> solid_;
  Probe probe_;

  void init();
  void apply_bc(VectorField2& vel);
  bool step(int n);

  void advect_diffuse();
  void build_divergence();
  bool pressure_poisson();
  void project();

  inline bool is_solid(int k) const { return solid_[k] != 0; }
};
