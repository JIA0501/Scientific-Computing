#include "FDSolver.hpp"
#include "Utils.hpp"
#include <cmath>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <omp.h>

FDSolver::FDSolver(const Config& cfg)
  : cfg_(cfg), g_(cfg.nx, cfg.ny, cfg.Lx, cfg.Ly),
    vel_(cfg.nx*cfg.ny), vel_tmp_(cfg.nx*cfg.ny),
    p_(cfg.nx*cfg.ny), div_(cfg.nx*cfg.ny)
{
  double D = 2.0*cfg_.R;
  cfg_.nu = cfg_.Uin * D / cfg_.Re;

  if(cfg_.dt <= 0){
    double h = std::min(g_.dx, g_.dy);
    cfg_.dt = 0.15 * h / std::max(1e-12, cfg_.Uin);
  }

  solid_ = utils::make_cylinder_mask(g_, cfg_.cx, cfg_.cy, cfg_.R);
}

void FDSolver::init(){
  #pragma omp parallel for
  for(int k=0;k<g_.nx*g_.ny;k++){
    vel_.u.a[k] = cfg_.Uin;
    vel_.v.a[k] = 0.0;
    p_.a[k] = 0.0;
  }

  #pragma omp parallel for
  for(int k=0;k<g_.nx*g_.ny;k++){
    if(is_solid(k)){
      vel_.u.a[k] = 0.0;
      vel_.v.a[k] = 0.0;
    }
  }
  apply_bc(vel_);

  double px = (cfg_.probe.x < 0) ? (cfg_.cx + 6.0*cfg_.R) : cfg_.probe.x;
  double py = (cfg_.probe.y < 0) ? cfg_.cy : cfg_.probe.y;
  std::ostringstream oss;
  oss << cfg_.out.dir << "/" << cfg_.probe.prefix << "_fd_Re" << (int)std::round(cfg_.Re) << ".csv";
  probe_.open(oss.str(), g_, px, py);
}

void FDSolver::apply_bc(VectorField2& vel){
  for(int i=0;i<g_.nx;i++){
    int kb = g_.idx(i,0);
    int kt = g_.idx(i,g_.ny-1);
    vel.u.a[kb]=0.0; vel.v.a[kb]=0.0;
    vel.u.a[kt]=0.0; vel.v.a[kt]=0.0;
  }

  for(int j=0;j<g_.ny;j++){
    int k = g_.idx(0,j);
    vel.u.a[k]=cfg_.Uin;
    vel.v.a[k]=0.0;
  }

  for(int j=0;j<g_.ny;j++){
    int kr = g_.idx(g_.nx-1,j);
    int km = g_.idx(g_.nx-2,j);
    vel.u.a[kr] = vel.u.a[km];
    vel.v.a[kr] = vel.v.a[km];
  }

  #pragma omp parallel for
  for(int k=0;k<g_.nx*g_.ny;k++){
    if(is_solid(k)){
      vel.u.a[k]=0.0; vel.v.a[k]=0.0;
    }
  }
}

void FDSolver::advect_diffuse(){
  const double dt = cfg_.dt, nu = cfg_.nu;
  const double dx = g_.dx, dy = g_.dy;
  const int nx = g_.nx, ny = g_.ny;

  #pragma omp parallel for collapse(2)
  for(int j=1;j<ny-1;j++){
    for(int i=1;i<nx-1;i++){
      int k = g_.idx(i,j);
      if(is_solid(k)){
        vel_tmp_.u.a[k]=0.0; vel_tmp_.v.a[k]=0.0;
        continue;
      }

      double u = vel_.u.a[k], v = vel_.v.a[k];

      double du_dx = (u > 0) ? (vel_.u.a[k] - vel_.u.a[g_.idx(i-1,j)]) / dx
                             : (vel_.u.a[g_.idx(i+1,j)] - vel_.u.a[k]) / dx;
      double du_dy = (v > 0) ? (vel_.u.a[k] - vel_.u.a[g_.idx(i,j-1)]) / dy
                             : (vel_.u.a[g_.idx(i,j+1)] - vel_.u.a[k]) / dy;

      double dv_dx = (u > 0) ? (vel_.v.a[k] - vel_.v.a[g_.idx(i-1,j)]) / dx
                             : (vel_.v.a[g_.idx(i+1,j)] - vel_.v.a[k]) / dx;
      double dv_dy = (v > 0) ? (vel_.v.a[k] - vel_.v.a[g_.idx(i,j-1)]) / dy
                             : (vel_.v.a[g_.idx(i,j+1)] - vel_.v.a[k]) / dy;

      double adv_u = u*du_dx + v*du_dy;
      double adv_v = u*dv_dx + v*dv_dy;

      double lap_u = (vel_.u.a[g_.idx(i+1,j)] - 2*vel_.u.a[k] + vel_.u.a[g_.idx(i-1,j)])/(dx*dx)
                   + (vel_.u.a[g_.idx(i,j+1)] - 2*vel_.u.a[k] + vel_.u.a[g_.idx(i,j-1)])/(dy*dy);
      double lap_v = (vel_.v.a[g_.idx(i+1,j)] - 2*vel_.v.a[k] + vel_.v.a[g_.idx(i-1,j)])/(dx*dx)
                   + (vel_.v.a[g_.idx(i,j+1)] - 2*vel_.v.a[k] + vel_.v.a[g_.idx(i,j-1)])/(dy*dy);

      vel_tmp_.u.a[k] = u + dt * (-adv_u + nu*lap_u);
      vel_tmp_.v.a[k] = v + dt * (-adv_v + nu*lap_v);
    }
  }

  apply_bc(vel_tmp_);
}

void FDSolver::build_divergence(){
  const double dx = g_.dx, dy = g_.dy;
  const int nx = g_.nx, ny = g_.ny;

  #pragma omp parallel for collapse(2)
  for(int j=1;j<ny-1;j++){
    for(int i=1;i<nx-1;i++){
      int k = g_.idx(i,j);
      if(is_solid(k)){ div_.a[k]=0.0; continue; }

      double du_dx = (vel_tmp_.u.a[g_.idx(i+1,j)] - vel_tmp_.u.a[g_.idx(i-1,j)])/(2.0*dx);
      double dv_dy = (vel_tmp_.v.a[g_.idx(i,j+1)] - vel_tmp_.v.a[g_.idx(i,j-1)])/(2.0*dy);
      div_.a[k] = du_dx + dv_dy;
    }
  }
}

bool FDSolver::pressure_poisson(){
  const double dx = g_.dx, dy = g_.dy;
  const double idx2 = 1.0/(dx*dx);
  const double idy2 = 1.0/(dy*dy);
  const double denom = 2.0*(idx2 + idy2);
  const double rhs_scale = cfg_.rho / cfg_.dt;
  const int nx = g_.nx, ny = g_.ny;

  // SOR 参数：1.5~1.9 通常都不错，先用 1.7
  const double omega = 1.7;

  bool converged = false;

  for(int it=0; it<cfg_.poissonIters; it++){
    double maxdiff = 0.0;

    // Gauss-Seidel: 原地更新 p_
    // 注意：这里不容易并行化（会破坏GS顺序），所以保持串行反而更稳
    for(int j=1;j<ny-1;j++){
      for(int i=1;i<nx-1;i++){
        int k = g_.idx(i,j);
        if(is_solid(k)){ p_.a[k]=0.0; continue; }

        double rhs = rhs_scale * div_.a[k];

        double pE = p_.a[g_.idx(i+1,j)];
        double pW = p_.a[g_.idx(i-1,j)];
        double pN = p_.a[g_.idx(i,j+1)];
        double pS = p_.a[g_.idx(i,j-1)];

        double p_gs = ( (pE+pW)*idx2 + (pN+pS)*idy2 - rhs ) / denom;

        double old = p_.a[k];
        double neu = (1.0-omega)*old + omega*p_gs;
        p_.a[k] = neu;

        double diff = std::abs(neu-old);
        if(diff > maxdiff) maxdiff = diff;
      }
    }

    // pressure BC:
    // inlet: dp/dx=0 ; outlet: p=0 ; walls: dp/dy=0
    for(int j=0;j<ny;j++){
      p_.a[g_.idx(0,j)] = p_.a[g_.idx(1,j)];
      p_.a[g_.idx(nx-1,j)] = 0.0;
    }
    for(int i=0;i<nx;i++){
      p_.a[g_.idx(i,0)] = p_.a[g_.idx(i,1)];
      p_.a[g_.idx(i,ny-1)] = p_.a[g_.idx(i,ny-2)];
    }

    if(maxdiff < cfg_.poissonTol){
      converged = true;
      break;
    }
  }

  return converged;
}

void FDSolver::project(){
  const double dt = cfg_.dt;
  const double dx = g_.dx, dy = g_.dy;
  const int nx=g_.nx, ny=g_.ny;

  #pragma omp parallel for collapse(2)
  for(int j=1;j<ny-1;j++){
    for(int i=1;i<nx-1;i++){
      int k=g_.idx(i,j);
      if(is_solid(k)){
        vel_.u.a[k]=0.0; vel_.v.a[k]=0.0;
        continue;
      }
      double dp_dx = (p_.a[g_.idx(i+1,j)] - p_.a[g_.idx(i-1,j)])/(2.0*dx);
      double dp_dy = (p_.a[g_.idx(i,j+1)] - p_.a[g_.idx(i,j-1)])/(2.0*dy);

      vel_.u.a[k] = vel_tmp_.u.a[k] - (dt/cfg_.rho) * dp_dx;
      vel_.v.a[k] = vel_tmp_.v.a[k] - (dt/cfg_.rho) * dp_dy;
    }
  }

  apply_bc(vel_);
}

bool FDSolver::step(int n){
  advect_diffuse();
  build_divergence();

  if(!pressure_poisson()){
  std::cerr << "[FD] WARN: Poisson not converged at step " << n << "\n";
  }
  project();

  if(utils::any_nonfinite(vel_)) return false;
  double ms = utils::max_speed(vel_);
  if(ms > cfg_.speed_cap_factor * cfg_.Uin) return false;

  if(cfg_.probe.every > 0 && (n % cfg_.probe.every)==0){
    double t = n * cfg_.dt;
    probe_.sample(t, g_, vel_);
  }

  if(cfg_.out.every > 0 && (n % cfg_.out.every)==0){
    std::ostringstream oss;
    oss << cfg_.out.dir << "/fd_step_" << n << ".csv";
    utils::write_csv_uvp_vort(oss.str(), g_, vel_, p_);
    std::cerr << "[FD] step " << n << " dt=" << cfg_.dt
              << " Re=" << cfg_.Re << " nu=" << cfg_.nu
              << " max|u|=" << ms << "\n";
  }

  return true;
}

bool FDSolver::run(){
  utils::ensure_dir(cfg_.out.dir);
  init();

  for(int n=0;n<cfg_.steps;n++){
    if(!step(n)){
      std::cerr << "[FD] BREAK at step " << n << " (Re=" << cfg_.Re << ")\n";
      probe_.close();
      return false;
    }
  }

  probe_.close();
  return true;
}
