#include "LBMSolver.hpp"
#include "Utils.hpp"
#include "Field.hpp"
#include <cmath>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <omp.h>

static constexpr int Q=9;
static constexpr int cx[Q]={0,1,0,-1,0,1,-1,-1,1};
static constexpr int cy[Q]={0,0,1,0,-1,1,1,-1,-1};
static constexpr double w[Q]={4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};

static inline double feq(int q,double rho,double ux,double uy){
  double cu = 3.0*(cx[q]*ux + cy[q]*uy);
  double uu = 1.5*(ux*ux + uy*uy);
  return w[q]*rho*(1.0 + cu + 0.5*cu*cu - uu);
}

LBMSolver::LBMSolver(const Config& cfg)
  : cfg_(cfg), g_(cfg.nx, cfg.ny, cfg.Lx, cfg.Ly)
{
  double scale_x = (cfg_.nx-1)/cfg_.Lx;
  double Dphys = 2.0*cfg_.R;
  double D_lbm = Dphys * scale_x;
  double nu_lbm = cfg_.U_lbm * D_lbm / cfg_.Re;

  cfg_.tau = (cfg_.tau > 0) ? cfg_.tau : (3.0*nu_lbm + 0.5);

  int N = g_.nx*g_.ny;
  f_.assign(N*Q, 0.0);
  f2_.assign(N*Q, 0.0);
  rho_.assign(N, 1.0);
  ux_.assign(N, 0.0);
  uy_.assign(N, 0.0);

  solid_ = utils::make_cylinder_mask(g_, cfg_.cx, cfg_.cy, cfg_.R);

  double px = (cfg_.probe.x < 0) ? (cfg_.cx + 6.0*cfg_.R) : cfg_.probe.x;
  double py = (cfg_.probe.y < 0) ? cfg_.cy : cfg_.probe.y;
  std::ostringstream oss;
  oss << cfg_.out.dir << "/" << cfg_.probe.prefix << "_lbm_Re" << (int)std::round(cfg_.Re) << ".csv";
  probe_.open(oss.str(), g_, px, py);
}

void LBMSolver::init(){
  #pragma omp parallel for
  for(int cell=0; cell<g_.nx*g_.ny; cell++){
    double rho=1.0;
    double ux=cfg_.U_lbm, uy=0.0;
    if(solid_[cell]){ ux=0.0; uy=0.0; }
    rho_[cell]=rho; ux_[cell]=ux; uy_[cell]=uy;
    int b=9*cell;
    for(int q=0;q<Q;q++) f_[b+q]=feq(q,rho,ux,uy);
  }
}

void LBMSolver::macroscopic(){
  #pragma omp parallel for
  for(int cell=0; cell<g_.nx*g_.ny; cell++){
    if(solid_[cell]){
      rho_[cell]=1.0; ux_[cell]=0.0; uy_[cell]=0.0;
      continue;
    }
    double rho=0.0, ux=0.0, uy=0.0;
    int b=9*cell;
    for(int q=0;q<Q;q++){
      double fq=f_[b+q];
      rho += fq;
      ux += fq*cx[q];
      uy += fq*cy[q];
    }
    ux /= rho;
    uy /= rho;
    rho_[cell]=rho;
    ux_[cell]=ux;
    uy_[cell]=uy;
  }
}

void LBMSolver::collide(){
  const double tau=cfg_.tau;
  const double om = 1.0/tau;

  #pragma omp parallel for
  for(int cell=0; cell<g_.nx*g_.ny; cell++){
    if(solid_[cell]) continue;
    int b=9*cell;
    double rho=rho_[cell], ux=ux_[cell], uy=uy_[cell];
    for(int q=0;q<Q;q++){
      double eq=feq(q,rho,ux,uy);
      f_[b+q] = f_[b+q] - om*(f_[b+q]-eq);
    }
  }
}

void LBMSolver::stream(){
  #pragma omp parallel for collapse(2)
  for(int j=0;j<g_.ny;j++){
    for(int i=0;i<g_.nx;i++){
      int cell=j*g_.nx+i;
      int b=9*cell;
      for(int q=0;q<Q;q++){
        int ip = i - cx[q];
        int jp = j - cy[q];
        if(ip<0 || ip>=g_.nx || jp<0 || jp>=g_.ny){
          f2_[b+q] = f_[b+q];
        } else {
          int from = 9*(jp*g_.nx+ip);
          f2_[b+q] = f_[from+q];
        }
      }
    }
  }

  auto bounce = [&](int cell){
    int b=9*cell;
    std::swap(f2_[b+1], f2_[b+3]);
    std::swap(f2_[b+2], f2_[b+4]);
    std::swap(f2_[b+5], f2_[b+7]);
    std::swap(f2_[b+6], f2_[b+8]);
  };

  #pragma omp parallel for
  for(int cell=0; cell<g_.nx*g_.ny; cell++){
    int j = cell / g_.nx;
    if(solid_[cell] || j==0 || j==g_.ny-1){
      bounce(cell);
    }
  }

  f_.swap(f2_);
}

void LBMSolver::apply_bc(){
  double U = cfg_.U_lbm;
  int nx=g_.nx;

  for(int j=1;j<g_.ny-1;j++){
    int cell = j*nx + 0;
    if(solid_[cell]) continue;
    int b=9*cell;

    double f0=f_[b+0], f2=f_[b+2], f4=f_[b+4], f3=f_[b+3], f6=f_[b+6], f7=f_[b+7];
    double rho = (f0+f2+f4 + 2.0*(f3+f6+f7)) / (1.0 - U);

    f_[b+1] = f3 + (2.0/3.0)*rho*U;
    f_[b+5] = f7 + 0.5*(f4-f2) + (1.0/6.0)*rho*U;
    f_[b+8] = f6 + 0.5*(f2-f4) + (1.0/6.0)*rho*U;
  }

  for(int j=1;j<g_.ny-1;j++){
    int c1=j*nx+(nx-1);
    int c0=j*nx+(nx-2);
    int b1=9*c1, b0=9*c0;
    for(int q=0;q<Q;q++) f_[b1+q]=f_[b0+q];
  }
}

void LBMSolver::write_out(int step){
  VectorField2 vel(g_.nx*g_.ny);
  ScalarField p(g_.nx*g_.ny);

  for(int k=0;k<g_.nx*g_.ny;k++){
    vel.u.a[k]=ux_[k];
    vel.v.a[k]=uy_[k];
    p.a[k]=rho_[k];
  }

  std::ostringstream oss;
  oss << cfg_.out.dir << "/lbm_step_" << step << ".csv";
  utils::write_csv_uvp_vort(oss.str(), g_, vel, p);

  double ms = utils::max_speed(vel);
  std::cerr << "[LBM] step " << step << " Re=" << cfg_.Re
            << " tau=" << cfg_.tau << " max|u|=" << ms << "\n";
}

bool LBMSolver::run(){
  utils::ensure_dir(cfg_.out.dir);
  init();

  for(int n=0;n<cfg_.steps;n++){
    macroscopic();
    collide();
    stream();
    apply_bc();

    if(n % 20 == 0){
      for(double r : rho_){
        if(!std::isfinite(r) || r <= 0.0){
          std::cerr << "[LBM] BREAK (bad rho) at step " << n << "\n";
          probe_.close();
          return false;
        }
      }
    }

    if(cfg_.probe.every > 0 && (n % cfg_.probe.every)==0){
      macroscopic();
      VectorField2 vel(g_.nx*g_.ny);
      for(int k=0;k<g_.nx*g_.ny;k++){
        vel.u.a[k]=ux_[k];
        vel.v.a[k]=uy_[k];
      }
      probe_.sample((double)n, g_, vel);

      double ms = utils::max_speed(vel);
      if(ms > cfg_.speed_cap_factor * cfg_.U_lbm){
        std::cerr << "[LBM] BREAK (speed cap) at step " << n << "\n";
        probe_.close();
        return false;
      }
    }

    if(cfg_.out.every > 0 && (n % cfg_.out.every)==0){
      macroscopic();
      write_out(n);
    }
  }

  probe_.close();
  return true;
}
