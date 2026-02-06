import pandas as pd
import matplotlib.pyplot as plt

# Read CSV (single column)
V = pd.read_csv("potential.csv", header=None)

# Generate x from index
x_V = V.index.values
y_V = V[0].values

fig, ax = plt.subplots(figsize=(10,5))

line_V, = ax.plot(x_V, y_V, label="Potential (V)")

ax.set_title("1D Poisson")
ax.set_xlabel("Grid Index")
ax.set_ylabel("Potential")
ax.legend()
ax.grid(alpha=0.3)

# --- Hover annotation ---
annot = ax.annotate(
    "",
    xy=(0,0),
    xytext=(30,30),
    textcoords="offset points",
    bbox=dict(boxstyle="round", fc="w"),
    arrowprops=dict(arrowstyle="->")
)
annot.set_visible(False)

def update_annot(line, ind):
    x, y = line.get_data()
    idx = ind["ind"][0]
    annot.xy = (x[idx], y[idx])
    annot.set_text(f"x={x[idx]}\ny={y[idx]:.4f}")

def hover(event):
    vis = annot.get_visible()
    if event.inaxes == ax:
        cont, ind = line_V.contains(event)
        if cont:
            update_annot(line_V, ind)
            annot.set_visible(True)
            fig.canvas.draw_idle()
            return
    if vis:
        annot.set_visible(False)
        fig.canvas.draw_idle()

fig.canvas.mpl_connect("motion_notify_event", hover)

plt.show()
