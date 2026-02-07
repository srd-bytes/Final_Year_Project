import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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



def plot_2d_potential_heatmap(csv_path, title="2D Poisson", cmap="viridis"):
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


def plot_2d_potential_3d(csv_path, title="2D Poisson Potential", cmap="viridis"):
    """
    Plot 2D potential as a 3D surface.
    CSV format: rows x cols (grid values)
    """

    V = pd.read_csv(csv_path, header=None).values
    Ny, Nx = V.shape

    x = np.arange(Nx)
    y = np.arange(Ny)
    X, Y = np.meshgrid(x, y)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")

    surf = ax.plot_surface(
        X, Y, V,
        cmap=cmap,
        linewidth=0,
        antialiased=True
    )

    ax.set_title(title)
    ax.set_xlabel("x index")
    ax.set_ylabel("y index")
    ax.set_zlabel("Potential (V)")

    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=12, label="Potential (V)")

    plt.show()

def plot_potential(csv_path, dim=1, **kwargs):
    if dim == 1:
        plot_1d_potential(csv_path, **kwargs)
    elif dim == 2:
        plot_2d_potential_3d(csv_path, **kwargs)
    else:
        raise ValueError("dim must be 1 or 2")


plot_potential("./result/potential_2D.csv", dim=2)
