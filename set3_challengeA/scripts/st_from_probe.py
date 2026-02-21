import numpy as np
import pandas as pd
import sys

if len(sys.argv) < 4:
    print("Usage: python st_from_probe.py <probe_csv> <Uin> <D>")
    sys.exit(1)

path = sys.argv[1]
Uin = float(sys.argv[2])
D = float(sys.argv[3])

df = pd.read_csv(path)
t = df["t"].to_numpy()
v = df["v"].to_numpy()

cut = int(0.30 * len(v))
t = t[cut:]
v = v[cut:] - np.mean(v[cut:])

dt = float(np.mean(np.diff(t)))
v = v * np.hanning(len(v))
V = np.fft.rfft(v)
freq = np.fft.rfftfreq(len(v), d=dt)

amp = np.abs(V)
mask = (freq >= 1.0) & (freq <= 5.0)  
amp2 = amp.copy()
amp2[~mask] = 0.0
k = int(np.argmax(amp2[1:]) + 1)

f = float(freq[k])
St = f * D / Uin

print(f"f={f:.6g}  St={St:.6g}  (dt_sample={dt:.6g}, N={len(v)})")
