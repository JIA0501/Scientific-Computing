# Assignment Set 3 – Challenge A (Kármán vortex street)

Implements the same cylinder-in-channel setup with:
1) Finite Difference (projection / fractional step) solver in C++ + OpenMP
2) Lattice Boltzmann (D2Q9 BGK) solver in C++ + OpenMP
3) Finite Element (ngsolve) solver in Python (fem/fem_ngsolve.py)

## Build
```bash
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -j
```

## Run (single Re)
FD:
```bash
./set3_chA --method fd --Re 100 --nx 800 --ny 200 --steps 20000 \
  --out dir=output,every=200 --probe every=10,prefix=probe
```

LBM:
```bash
./set3_chA --method lbm --Re 200 --nx 1200 --ny 300 --steps 40000 \
  --out dir=output,every=400 --probe every=20,prefix=probe --Ulbm 0.08
```

FEM (requires ngsolve installed):
```bash
./set3_chA --method fem --Re 100 --steps 2000
```

## Sweep Re (auto-find max stable Re)
```bash
./set3_chA --method fd  --sweep start=50,growth=1.25,tol=10,steps=8000
./set3_chA --method lbm --sweep start=50,growth=1.25,tol=10,steps=15000
```

## Strouhal number from probe
Probe output is written to `output/probe_<method>_ReXXX.csv`.
Compute St:
```bash
python3 scripts/st_from_probe.py output/probe_fd_Re100.csv 1.0 0.1
```
Where `Uin=1.0`, `D=2R=0.1` for the default geometry.
