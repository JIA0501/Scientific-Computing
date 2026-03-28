#pragma once
#include <string>
#include <fstream>
#include "Grid.hpp"
#include "Field.hpp"

class Probe {
public:
  Probe() = default;
  void open(const std::string& path, const Grid& g, double x, double y);
  void sample(double t, const Grid& g, const VectorField2& vel);
  void close();

  int i() const { return pi_; }
  int j() const { return pj_; }

private:
  int pi_ = -1, pj_ = -1;
  std::ofstream f_;
};
