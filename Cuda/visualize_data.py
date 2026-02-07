import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

def plot_1d_potential(csv_path, title="1D Poisson"):
    """
    Plot 1D potential from CSV with hover annotation.
    CSV format: single column of values.
    """

    V = pd.read_csv(csv_path, header=None)

    x = V.index.values
    y = V[0].values

    fig, ax = plt.subplots(figsize=(10, 5))
    line, = ax.plot(x, y, label="Potential (V)")

    ax.set_title(title)
    ax.set_xlabel("Grid Index")
    ax.set_ylabel("Potential")
    ax.legend()
    ax.grid(alpha=0.3)

    # Hover annotation
    annot = ax.annotate(
        "",
        xy=(0, 0),
        xytext=(30, 30),
        textcoords="offset points",
        bbox=dict(boxstyle="round", fc="w"),
        arrowprops=dict(arrowstyle="->")
    )
    annot.set_visible(False)

    def update_annot(ind):
        idx = ind["ind"][0]
        annot.xy = (x[idx], y[idx])
        annot.set_text(f"x={x[idx]}\nV={y[idx]:.6f}")

    def hover(event):
        if event.inaxes == ax:
            cont, ind = line.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
                return
        annot.set_visible(False)
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", hover)
    plt.show()



def plot_2d_potential_heatmap(csv_path, title="2D Poisson", cmap="rainbow"):
    """
    Plot 2D potential as a heatmap.
    CSV format: 2D grid (rows x cols)
    """

    V = pd.read_csv(csv_path, header=None).values

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(V, origin="lower", cmap=cmap, aspect="auto")

    ax.set_title(title)
    ax.set_xlabel("x index")
    ax.set_ylabel("y index")

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Potential (V)")

    plt.show()







def animate_3d_heatmap_slices(
    csv_path,
    interval=5,
    cmap="viridis"
):
    """
    Animate Z-slices of a 3D scalar field as a 2D heatmap.
    """
    slices = []
    current_slice = []

    with open(csv_path, "r") as f:
        for line in f:
            line = line.strip()

            # New z-slice marker
            if line.startswith("#"):
                if current_slice:
                    slices.append(np.array(current_slice, dtype=float))
                    current_slice = []
                continue

            # Skip empty lines
            if not line:
                continue

            # Data row
            row = list(map(float, line.split(",")))
            current_slice.append(row)

        # Append last slice
        if current_slice:
            slices.append(np.array(current_slice, dtype=float))

    phi = np.stack(slices, axis=0)  # (Nz, Ny, Nx)
    Nz, Ny, Nx = phi.shape

    fig, ax = plt.subplots(figsize=(6, 5))

    # Global color scale (IMPORTANT)
    vmin = phi.min()
    vmax = phi.max()

    im = ax.imshow(
        phi[0],
        origin="lower",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        animated=True
    )

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Potential (V)")

    title = ax.set_title("z = 0")

    def update(z):
        im.set_array(phi[z])
        title.set_text(f"z = {z}")
        return im, title

    anim = FuncAnimation(
        fig,
        update,
        frames=Nz,
        interval=interval,
        blit=False
    )

    ax.set_xlabel("X index")
    ax.set_ylabel("Y index")
    ax.set_aspect("equal")

    plt.show()


# --------------------------Test---------------------------------------

# plot_1d_potential("./result/potential_1D.csv")
# plot_2d_potential_heatmap("./result/potential_2D_new.csv")
animate_3d_heatmap_slices("./result/potential_3D.csv", interval=2, cmap="rainbow")


