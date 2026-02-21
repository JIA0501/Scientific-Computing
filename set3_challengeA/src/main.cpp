#include "Config.hpp"
#include "FDSolver.hpp"
#include "LBMSolver.hpp"
#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <functional>
#include <omp.h>

static void usage(){
  std::cerr <<
  "Challenge A: Karman vortex street (FD / LBM / FEM)\n"
  "Usage:\n"
  "  ./set3_chA --method fd|lbm|fem --Re <val> [--nx N --ny N --steps S]\n"
  "Sweep:\n"
  "  ./set3_chA --method fd|lbm --sweep start=<Re0>,growth=<g>,tol=<tol>,steps=<steps>\n"
  "Output:\n"
  "  --out dir=<path>,every=<k>\n"
  "Probe:\n"
  "  --probe every=<k>,x=<val>,y=<val>,prefix=<name>\n";
}

static bool parse_kvlist(const std::string& s,
                         const std::function<void(const std::string&, const std::string&)>& fn){
  size_t start=0;
  while(start < s.size()){
    size_t end = s.find(',', start);
    if(end==std::string::npos) end = s.size();
    std::string item = s.substr(start, end-start);
    size_t eq = item.find('=');
    if(eq==std::string::npos) return false;
    std::string k = item.substr(0,eq);
    std::string v = item.substr(eq+1);
    fn(k,v);
    start = end + 1;
  }
  return true;
}

static bool run_once(Config cfg){
  if(cfg.threads > 0) omp_set_num_threads(cfg.threads);

  if(cfg.method=="fd"){
    FDSolver s(cfg);
    return s.run();
  } else if(cfg.method=="lbm"){
    LBMSolver s(cfg);
    return s.run();
  } else if(cfg.method=="fem"){
    std::string cmd = "python3 fem/fem_ngsolve.py --Re " + std::to_string(cfg.Re)
                    + " --steps " + std::to_string(cfg.steps);
    std::cerr << "[FEM] running: " << cmd << "\n";
    int rc = std::system(cmd.c_str());
    return rc==0;
  } else {
    std::cerr << "Unknown method: " << cfg.method << "\n";
    return false;
  }
}

static int do_sweep(Config cfg){
  if(cfg.method=="fem"){
    std::cerr << "Sweep for FEM is not enabled in this template.\n";
    return 1;
  }

  double Re = cfg.sweep.startRe;
  double last_ok = Re;
  double first_bad = -1.0;
  int runs = 0;

  auto run_at = [&](double ReVal)->bool{
    Config c = cfg;
    c.Re = ReVal;
    c.steps = cfg.sweep.warmup_steps;
    std::cerr << "\n=== SWEEP RUN: method=" << c.method << " Re=" << c.Re
              << " steps=" << c.steps << " ===\n";
    return run_once(c);
  };

  while(runs < cfg.sweep.maxRuns){
    runs++;
    bool ok = run_at(Re);
    if(ok){
      last_ok = Re;
      Re *= cfg.sweep.growth;
    } else {
      first_bad = Re;
      break;
    }
  }

  if(first_bad < 0){
    std::cerr << "Sweep reached maxRuns without failure. last_ok=" << last_ok << "\n";
    return 0;
  }

  while((first_bad - last_ok) > cfg.sweep.tol && runs < cfg.sweep.maxRuns){
    runs++;
    double mid = 0.5*(first_bad + last_ok);
    bool ok = run_at(mid);
    if(ok) last_ok = mid;
    else   first_bad = mid;
  }

  std::cerr << "\n=== SWEEP RESULT ===\n";
  std::cerr << "method=" << cfg.method
            << " max_stable_Re≈" << last_ok
            << " (bracket [" << last_ok << ", " << first_bad << "], tol=" << cfg.sweep.tol << ")\n";
  return 0;
}

int main(int argc, char** argv){
  Config cfg;

  auto need = [&](int& i, const char* name)->std::string{
    if(i+1>=argc){ std::cerr << "Missing value for " << name << "\n"; usage(); std::exit(1); }
    return std::string(argv[++i]);
  };

  for(int i=1;i<argc;i++){
    std::string a=argv[i];

    if(a=="--help" || a=="-h"){ usage(); return 0; }
    else if(a=="--method"){ cfg.method = need(i,"--method"); }
    else if(a=="--Re"){ cfg.Re = std::stod(need(i,"--Re")); }
    else if(a=="--nx"){ cfg.nx = std::stoi(need(i,"--nx")); }
    else if(a=="--ny"){ cfg.ny = std::stoi(need(i,"--ny")); }
    else if(a=="--steps"){ cfg.steps = std::stoi(need(i,"--steps")); }
    else if(a=="--Uin"){ cfg.Uin = std::stod(need(i,"--Uin")); }
    else if(a=="--dt"){ cfg.dt = std::stod(need(i,"--dt")); }
    else if(a=="--threads"){ cfg.threads = std::stoi(need(i,"--threads")); }
    else if(a=="--poissonIters"){ cfg.poissonIters = std::stoi(need(i,"--poissonIters")); }
    else if(a=="--poissonTol"){ cfg.poissonTol = std::stod(need(i,"--poissonTol")); }
    else if(a=="--speedcap"){ cfg.speed_cap_factor = std::stod(need(i,"--speedcap")); }
    else if(a=="--Ulbm"){ cfg.U_lbm = std::stod(need(i,"--Ulbm")); }
    else if(a=="--tau"){ cfg.tau = std::stod(need(i,"--tau")); }
    else if(a=="--out"){
      std::string v = need(i,"--out");
      bool ok = parse_kvlist(v, [&](const std::string& k, const std::string& val){
        if(k=="dir") cfg.out.dir = val;
        else if(k=="every") cfg.out.every = std::stoi(val);
      });
      if(!ok){ std::cerr << "Bad --out format\n"; usage(); return 1; }
    }
    else if(a=="--probe"){
      std::string v = need(i,"--probe");
      bool ok = parse_kvlist(v, [&](const std::string& k, const std::string& val){
        if(k=="every") cfg.probe.every = std::stoi(val);
        else if(k=="x") cfg.probe.x = std::stod(val);
        else if(k=="y") cfg.probe.y = std::stod(val);
        else if(k=="prefix") cfg.probe.prefix = val;
      });
      if(!ok){ std::cerr << "Bad --probe format\n"; usage(); return 1; }
    }
    else if(a=="--sweep"){
      cfg.sweep.enabled = true;
      std::string v = need(i,"--sweep");
      bool ok = parse_kvlist(v, [&](const std::string& k, const std::string& val){
        if(k=="start") cfg.sweep.startRe = std::stod(val);
        else if(k=="growth") cfg.sweep.growth = std::stod(val);
        else if(k=="tol") cfg.sweep.tol = std::stod(val);
        else if(k=="steps") cfg.sweep.warmup_steps = std::stoi(val);
        else if(k=="maxruns") cfg.sweep.maxRuns = std::stoi(val);
      });
      if(!ok){ std::cerr << "Bad --sweep format\n"; usage(); return 1; }
    }
    else {
      std::cerr << "Unknown arg: " << a << "\n";
      usage();
      return 1;
    }
  }

  if(cfg.sweep.enabled) return do_sweep(cfg);

  bool ok = run_once(cfg);
  return ok ? 0 : 2;
}
