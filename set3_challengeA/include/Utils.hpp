#pragma once
#include <string>
#include <vector>
#include <cstdint>
#include "Grid.hpp"
#include "Field.hpp"

namespace utils {
  void ensure_dir(const std::string& dir);

  void write_csv_uvp_vort(const std::string& path,
                          const Grid& g,
                          const VectorField2& vel,
                          const ScalarField& p);

  ScalarField compute_vorticity(const Grid& g, const VectorField2& vel);

  std::vector<uint8_t> make_cylinder_mask(const Grid& g, double cx, double cy, double R);

  bool any_nonfinite(const VectorField2& vel);
  double max_speed(const VectorField2& vel);
}
