import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ====== parameter ======
N = 100                
Nx, Ny = N, N + 1
times = [0.0, 0.001, 0.01, 0.1, 1.0]

BASE = Path(__file__).parent


OUT = BASE / "plots"
OUT.mkdir(exist_ok=True)

# ------------------------------------------------------------
# E: profile plots (y, numerical, analytic)
# ------------------------------------------------------------
def plot_E_profiles():
    for t in times:
        fn = BASE / f"profile_t_{t:.6f}.dat"
        if not fn.exists():
            print(f"[E] skip (not found): {fn.name}")
            continue

        data = np.loadtxt(fn)
        y, num, ana = data[:, 0], data[:, 1], data[:, 2]

        plt.figure()
        plt.plot(y, num, label="numerical")
        plt.plot(y, ana, label="analytic")
        plt.xlabel("y")
        plt.ylabel("c(y,t)")
        plt.title(f"Profile at t={t}")
        plt.legend()
        plt.tight_layout()

        out_png = OUT / f"E_profile_t_{t:.6f}.png"
        plt.savefig(out_png, dpi=200)
        plt.show()
        plt.close()
        print(f"[E] saved: {out_png.name}")


# ------------------------------------------------------------
# F: 2D field heatmaps (x, y, c)
# ------------------------------------------------------------
def plot_F_fields():
    for t in times:
        fn = BASE / f"field_t_{t:.6f}.xyz"
        if not fn.exists():
            print(f"[F] skip (not found): {fn.name}")
            continue

        data = np.loadtxt(fn)
        # data columns: x, y, c
        c = data[:, 2]

        # -> reshape 成 (Ny, Nx)
        C = c.reshape(Ny, Nx)

        plt.figure()
        im = plt.imshow(
            C, origin="lower",
            extent=[0, 1, 0, 1],
            aspect="auto"
        )
        plt.colorbar(im, label="c")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title(f"2D concentration field at t={t}")
        plt.tight_layout()

        out_png = OUT / f"F_field_t_{t:.6f}.png"
        plt.savefig(out_png, dpi=200)
        plt.show()
        plt.close()
        print(f"[F] saved: {out_png.name}")


# ------------------------------------------------------------
# G: animated plot until equilibrium
# ------------------------------------------------------------
def make_G_animation_gif():
    try:
        import matplotlib.animation as animation
    except Exception as e:
        print("matplotlib.animation import failed:", e)
        return

    
    field_files = sorted(BASE.glob("field_t_*.xyz"))
    if len(field_files) < 2:
        print("[G] need at least 2 field_t_*.xyz files to animate.")
        return

    data0 = np.loadtxt(field_files[0])
    C0 = data0[:, 2].reshape(Ny, Nx)

    fig, ax = plt.subplots()
    im = ax.imshow(C0, origin="lower", extent=[0,1,0,1], aspect="auto")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("c")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    title = ax.set_title(field_files[0].stem)

    def update(frame_idx):
        fn = field_files[frame_idx]
        data = np.loadtxt(fn)
        C = data[:, 2].reshape(Ny, Nx)
        im.set_data(C)
        title.set_text(fn.stem)
        return (im, title)

    ani = animation.FuncAnimation(fig, update, frames=len(field_files), interval=400, blit=False)

    out_gif = OUT / "G_diffusion_animation.gif"
    ani.save(out_gif, writer="pillow", dpi=150)
    plt.close(fig)
    print(f"[G] saved: {out_gif.name}")


if __name__ == "__main__":
    plot_E_profiles()
    plot_F_fields()
    make_G_animation_gif()
