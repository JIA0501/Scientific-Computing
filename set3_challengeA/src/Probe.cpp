#include "Probe.hpp"
#include <cmath>
#include <filesystem>

void Probe::open(const std::string& path, const Grid& g, double x, double y){
  std::filesystem::create_directories(std::filesystem::path(path).parent_path());
  pi_ = (int)std::round(x / g.dx);
  pj_ = (int)std::round(y / g.dy);
  if(pi_ < 0) pi_ = 0;
  if(pi_ >= g.nx) pi_ = g.nx-1;
  if(pj_ < 0) pj_ = 0;
  if(pj_ >= g.ny) pj_ = g.ny-1;

  f_.open(path);
  f_ << "t,u,v,i,j\n";
}

void Probe::sample(double t, const Grid& g, const VectorField2& vel){
  if(!f_.is_open()) return;
  int k = g.idx(pi_, pj_);
  f_ << t << "," << vel.u.a[k] << "," << vel.v.a[k] << "," << pi_ << "," << pj_ << "\n";
}

void Probe::close(){
  if(f_.is_open()) f_.close();
}
