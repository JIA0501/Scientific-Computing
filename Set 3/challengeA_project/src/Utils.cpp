#include "Utils.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace fs = std::filesystem;

void utils::ensure_dir(const std::string& dir){
  fs::create_directories(dir);
}

ScalarField utils::compute_vorticity(const Grid& g, const VectorField2& vel){
  ScalarField w(g.nx*g.ny);
  for(int j=1;j<g.ny-1;j++){
    for(int i=1;i<g.nx-1;i++){
      int k=g.idx(i,j);
      double dv_dx = (vel.v[g.idx(i+1,j)] - vel.v[g.idx(i-1,j)])/(2.0*g.dx);
      double du_dy = (vel.u[g.idx(i,j+1)] - vel.u[g.idx(i,j-1)])/(2.0*g.dy);
      w[k]=dv_dx - du_dy;
    }
  }
  return w;
}

void utils::write_csv_uvp_vort(const std::string& path, const Grid& g,
                               const VectorField2& vel, const ScalarField& p){
  ScalarField w = compute_vorticity(g, vel);
  std::ofstream f(path);
  f << "i,j,x,y,u,v,p,vort\n";
  f << std::setprecision(10);
  for(int j=0;j<g.ny;j++){
    for(int i=0;i<g.nx;i++){
      int k=g.idx(i,j);
      double x=i*g.dx, y=j*g.dy;
      f << i << "," << j << "," << x << "," << y << ","
        << vel.u.a[k] << "," << vel.v.a[k] << "," << p.a[k] << "," << w.a[k] << "\n";
    }
  }
}

std::vector<uint8_t> utils::make_cylinder_mask(const Grid& g, double cx, double cy, double R){
  std::vector<uint8_t> m(g.nx*g.ny, 0);
  double R2 = R*R;
  for(int j=0;j<g.ny;j++){
    for(int i=0;i<g.nx;i++){
      double x=i*g.dx, y=j*g.dy;
      double d2 = (x-cx)*(x-cx) + (y-cy)*(y-cy);
      if(d2 <= R2) m[g.idx(i,j)] = 1;
    }
  }
  return m;
}

bool utils::any_nonfinite(const VectorField2& vel){
  for(int k=0;k<vel.u.n;k++){
    if(!std::isfinite(vel.u.a[k]) || !std::isfinite(vel.v.a[k])) return true;
  }
  return false;
}

double utils::max_speed(const VectorField2& vel){
  double m = 0.0;
  for(int k=0;k<vel.u.n;k++){
    double s = std::sqrt(vel.u.a[k]*vel.u.a[k] + vel.v.a[k]*vel.v.a[k]);
    if(s > m) m = s;
  }
  return m;
}
