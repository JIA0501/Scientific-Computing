import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ====== choose one file to animate ======
# filename = r"A:\science_computation\A1_1\wave_i.txt"
# filename = r"A:\science_computation\A1_1\wave_ii.txt"
filename = r"A:\science_computation\A1_1\wave_iii.txt"


def read_blocks(filename):
    """
    Read blocks:
      # t = ...
      x psi
    Return: times(list[float]), X(list[list[float]]), PSI(list[list[float]])
    """
    times = []
    X = []
    PSI = []
    cur_x, cur_psi = [], []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("#"):
                # store previous block
                if cur_x:
                    X.append(cur_x)
                    PSI.append(cur_psi)
                    cur_x, cur_psi = [], []
                # parse time number after '='
                # e.g. "# t = 0.501"
                t = float(line.split("=")[1])
                times.append(t)
            else:
                x, psi = map(float, line.split())
                cur_x.append(x)
                cur_psi.append(psi)

    if cur_x:
        X.append(cur_x)
        PSI.append(cur_psi)

    return times, X, PSI


times, X, PSI = read_blocks(filename)

# ====== set up figure ======
fig, ax = plt.subplots()
line, = ax.plot(X[0], PSI[0])
ax.set_xlabel("x")
ax.set_ylabel("psi")

# fix y-limits based on global data range for nicer animation
all_vals = [v for block in PSI for v in block]
ymin, ymax = min(all_vals), max(all_vals)
pad = 0.05 * (ymax - ymin + 1e-12)
ax.set_ylim(ymin - pad, ymax + pad)

title = ax.set_title(f"t = {times[0]:.3f}")


def update(frame):
    line.set_data(X[frame], PSI[frame])
    title.set_text(f"t = {times[frame]:.3f}")
    return line, title


ani = FuncAnimation(fig, update, frames=len(times), interval=60, blit=True)

plt.show()

# ====== (optional) save animation ======
# GIF (needs pillow): pip install pillow
# ani.save("wave_animation.gif", writer="pillow", fps=20)

# MP4 (needs ffmpeg installed and in PATH)
# ani.save("wave_animation.mp4", writer="ffmpeg", fps=20)
