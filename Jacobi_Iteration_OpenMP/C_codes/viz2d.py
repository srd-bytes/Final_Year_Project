import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --- Load CSV containing only V column ---
df = pd.read_csv("potential_2d.csv", header=None, names=["V"])

# Infer grid size n × n
N = int(np.sqrt(len(df)))

# Convert linear 1D array → 2D grid
# First N values correspond to i=0, j=0..N-1
V = df["V"].values.reshape(N, N)

# Create mesh for plotting
x = np.arange(N)
y = np.arange(N)
X, Y = np.meshgrid(x, y)

# --- 3D Surface Plot ---
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(
    X, Y, V,
    cmap='plasma',
    edgecolor='none'
)

fig.colorbar(surf, shrink=0.6, aspect=10, label="Potential V")

ax.set_title("3D Surface Plot of 2D Poisson Potential")
ax.set_xlabel("j index (column)")
ax.set_ylabel("i index (row)")
ax.set_zlabel("Potential V")

plt.show()
