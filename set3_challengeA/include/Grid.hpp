#pragma once

struct Grid {
  int nx, ny;
  double Lx, Ly;
  double dx, dy;

  Grid(int nx_, int ny_, double Lx_, double Ly_)
    : nx(nx_), ny(ny_), Lx(Lx_), Ly(Ly_) {
    dx = Lx / (nx - 1);
    dy = Ly / (ny - 1);
  }

  inline int idx(int i, int j) const { return j * nx + i; }
  inline bool inside(int i, int j) const { return (i>=0 && i<nx && j>=0 && j<ny); }
};
