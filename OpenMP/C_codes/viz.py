import pandas as pd
import matplotlib.pyplot as plt

rho = pd.read_csv("charge_1d.csv", header=None)
V = pd.read_csv("potential_1d.csv", header=None)

x_rho, y_rho = rho[0], rho[1]
x_V, y_V = V[0], V[1]

fig, ax = plt.subplots(figsize=(10,5))

line_rho, = ax.plot(x_rho, y_rho, label="Charge (ρ)")
line_V,   = ax.plot(x_V, y_V, label="Potential (V)")

ax.set_title("1D Poisson")
ax.set_xlabel("Grid Index")
ax.legend()
ax.grid(alpha=0.3)

# --- Hover annotation ---
annot = ax.annotate("", xy=(0,0), xytext=(30,30),
                    textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->"))
annot.set_visible(False)

def update_annot(line, ind):
    x, y = line.get_data()
    index = ind["ind"][0]
    annot.xy = (x[index], y[index])
    text = f"x={x[index]:.2f}\ny={y[index]:.2f}"
    annot.set_text(text)

def hover(event):
    vis = annot.get_visible()
    if event.inaxes == ax:
        for line in [line_rho, line_V]:
            cont, ind = line.contains(event)
            if cont:
                update_annot(line, ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
                return
    if vis:
        annot.set_visible(False)
        fig.canvas.draw_idle()

fig.canvas.mpl_connect("motion_notify_event", hover)

plt.show()
