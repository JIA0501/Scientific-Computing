import os, re, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

Uin = 1.0
D = 0.1

# --- same logic as st_from_probe.py ---
def st_from_probe_csv(path, Uin=1.0, D=0.1, cut_frac=0.30):
    df = pd.read_csv(path)
    if "t" not in df.columns or "v" not in df.columns:
        raise ValueError(f"Missing columns in {path}. Need t and v.")
    t = df["t"].to_numpy(dtype=float)
    v = df["v"].to_numpy(dtype=float)

    if len(v) < 16:
        return np.nan, np.nan, np.nan, len(v)

    cut = int(cut_frac * len(v))
    t = t[cut:]
    v = v[cut:]
    v = v - np.mean(v)

    dt = float(np.mean(np.diff(t)))
    if not np.isfinite(dt) or dt <= 0:
        return np.nan, np.nan, np.nan, len(v)

    V = np.fft.rfft(v)
    freq = np.fft.rfftfreq(len(v), d=dt)

    # ignore DC, take max peak
    k = int(np.argmax(np.abs(V[1:])) + 1)
    f = float(freq[k])
    St = f * D / Uin
    return f, St, dt, len(v)

def extract_re(filename):
    m = re.search(r"Re(\d+)", os.path.basename(filename))
    return int(m.group(1)) if m else None

def collect(method, pattern):
    rows = []
    for fn in glob.glob(pattern):
        Re = extract_re(fn)
        if Re is None:
            continue
        f, St, dt, N = st_from_probe_csv(fn, Uin=Uin, D=D)
        rows.append({"method": method, "Re": Re, "f": f, "St": St, "dt": dt, "N": N, "file": fn})
    rows.sort(key=lambda r: r["Re"])
    return rows

def main():
    # Run from build folder typically. We support relative paths.
    # FD probes (you already have: build/output/probe_fd_ReXXX.csv)
    fd_rows  = collect("FD",  os.path.join("output", "probe_fd_Re*.csv"))

    # FEM probes: build/../fem/fem_output/probe_fem_ReXXX.csv
    fem_rows = collect("FEM", os.path.join("..", "fem", "fem_output", "probe_fem_Re*.csv"))

    # LBM probes: build/output_lbm/probe_lbm_ReXXX.csv
    lbm_rows = collect("LBM", os.path.join("output_lbm", "probe_lbm_Re*.csv"))

    all_rows = fd_rows + fem_rows + lbm_rows
    if not all_rows:
        raise SystemExit("No probe files found. Check your directories/patterns.")

    df = pd.DataFrame(all_rows)
    out_csv = "st_summary.csv"
    df.to_csv(out_csv, index=False)
    print(f"[OK] Wrote {out_csv} with {len(df)} rows")

    # Plot
    plt.figure(figsize=(7.5, 4.8))

    for method, marker in [("FD", "o"), ("FEM", "s"), ("LBM", "^")]:
        sub = df[df["method"] == method].sort_values("Re")
        if len(sub) == 0:
            continue
        plt.plot(sub["Re"], sub["St"], marker=marker, linewidth=1.5, label=method)

    plt.xscale("log")
    plt.yscale("log")  # important for your current results
    plt.xlabel("Reynolds number Re")
    plt.ylabel("Strouhal number St = f D / U")
    plt.title("Cylinder wake: Strouhal number vs Reynolds number (FD vs FEM vs LBM)")
    plt.grid(True, which="both", linestyle="--", linewidth=0.6)
    plt.legend()
    plt.tight_layout()

    out_png = "st_vs_re.png"
    plt.savefig(out_png, dpi=300)
    print(f"[OK] Wrote {out_png}")

if __name__ == "__main__":
    main()

