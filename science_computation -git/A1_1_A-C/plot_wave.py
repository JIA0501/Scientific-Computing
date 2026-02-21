import matplotlib.pyplot as plt


filename1 = r"A:\science_computation\A1_1\wave_i.txt"
filename2 = r"A:\science_computation\A1_1\wave_ii.txt"
filename3 = r"A:\science_computation\A1_1\wave_iii.txt"


def read_wave_file(filename):
    """
    Read the C++ output format:
      # t = ...
      x psi
      x psi
      ...
    (blank lines separate time blocks)

    Returns:
      times (list[str]),
      x_vals (list[list[float]]),
      psi_vals (list[list[float]])
    """
    times = []
    x_vals = []
    psi_vals = []

    cur_x = []
    cur_psi = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("#"):
                # New time block: store the previous block (if any)
                if cur_x:
                    x_vals.append(cur_x)
                    psi_vals.append(cur_psi)
                    cur_x = []
                    cur_psi = []
                times.append(line.replace("#", "").strip())  # e.g., "t = 0.501"
            else:
                x, psi = map(float, line.split())
                cur_x.append(x)
                cur_psi.append(psi)

    # Finalize: store the last block
    if cur_x:
        x_vals.append(cur_x)
        psi_vals.append(cur_psi)

    return times, x_vals, psi_vals


def plot_case(filename, title):
    times, x_vals, psi_vals = read_wave_file(filename)

    plt.figure()
    for i in range(len(x_vals)):
        # In case the number of time labels doesn't match the number of blocks
        label = times[i] if i < len(times) else f"snapshot {i}"
        plt.plot(x_vals[i], psi_vals[i], label=label)

    plt.xlabel("x")
    plt.ylabel("psi")
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


# ====== Plot the three initial-condi
plot_case(filename1, "Wave evolution: i) ψ(x,0)=sin(2πx)")
plot_case(filename2, "Wave evolution: ii) ψ(x,0)=sin(5πx)")
plot_case(filename3, "Wave evolution: iii) ψ(x,0)=sin(5πx) for 1/5<x<2/5 else 0")