#pragma once
#include <vector>

struct ScalarField {
  int n = 0;
  std::vector<double> a;
  ScalarField() = default;
  explicit ScalarField(int n_) : n(n_), a(n_, 0.0) {}
  inline double& operator[](int k){ return a[k]; }
  inline const double& operator[](int k) const { return a[k]; }
};

struct VectorField2 {
  ScalarField u, v;
  VectorField2() = default;
  explicit VectorField2(int n) : u(n), v(n) {}
};
