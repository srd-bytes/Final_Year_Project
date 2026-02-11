import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.style.use("seaborn-v0_8-darkgrid")

# -------------------------------------------------------
# Helper Function
# -------------------------------------------------------
def load_and_merge(cpu_file, gpu_file):

    cpu = pd.read_csv(cpu_file)
    gpu = pd.read_csv(gpu_file)

    data = pd.merge(cpu, gpu, on="Power", how="outer")
    data = data.sort_values("Power")

    data["Speedup"] = data["CPU"] / data["GPU"]
    data["N"] = 2 ** data["Power"]

    return data

def annotate_points(ax, x, y, labels, offset=(5,5)):
    for xi, yi, lab in zip(x, y, labels):
        if not np.isnan(yi):
            ax.annotate(f"{int(lab)}",
                        (xi, yi),
                        textcoords="offset points",
                        xytext=offset,
                        fontsize=8)

# -------------------------------------------------------
# Load All Dimensions
# -------------------------------------------------------
data_1d = load_and_merge("./Cpu/result/performance/benchmark_cpu_1d.csv", "./Cuda/result/performance/benchmark_gpu_1d.csv")
data_2d = load_and_merge("./Cpu/result/performance/benchmark_cpu_2d.csv", "./Cuda/result/performance/benchmark_gpu_2d.csv")
data_3d = load_and_merge("./Cpu/result/performance/benchmark_cpu_3d.csv", "./Cuda/result/performance/benchmark_gpu_3d.csv")


# -------------------------------------------------------
# 1️⃣ CPU vs GPU Scaling (All Dimensions)
# -------------------------------------------------------
fig, ax = plt.subplots(figsize=(10,7))

ax.loglog(data_1d["N"], data_1d["GPU"], 'o-', label="GPU 1D")
ax.loglog(data_2d["N"], data_2d["GPU"], 's-', label="GPU 2D")
ax.loglog(data_3d["N"], data_3d["GPU"], '^-', label="GPU 3D")

ax.loglog(data_1d["N"], data_1d["CPU"], '--', alpha=0.6)
ax.loglog(data_2d["N"], data_2d["CPU"], '--', alpha=0.6)
ax.loglog(data_3d["N"], data_3d["CPU"], '--', alpha=0.6)

# Annotate GPU points
annotate_points(ax, data_1d["N"], data_1d["GPU"], data_1d["Power"])
annotate_points(ax, data_2d["N"], data_2d["GPU"], data_2d["Power"])
annotate_points(ax, data_3d["N"], data_3d["GPU"], data_3d["Power"])

ax.set_xlabel("Grid Size N")
ax.set_ylabel("Time (seconds)")
ax.set_title("Jacobi Scaling: CPU vs GPU (1D / 2D / 3D)")
ax.legend()

plt.tight_layout()
plt.savefig("scaling_all_dimensions.png", dpi=600)
plt.show()


# -------------------------------------------------------
# 2️⃣ Speedup Comparison
# -------------------------------------------------------
fig, ax = plt.subplots(figsize=(10,7))

ax.semilogx(data_1d["N"], data_1d["Speedup"], 'o-', label="1D")
ax.semilogx(data_2d["N"], data_2d["Speedup"], 's-', label="2D")
ax.semilogx(data_3d["N"], data_3d["Speedup"], '^-', label="3D")

annotate_points(ax, data_1d["N"], data_1d["Speedup"], data_1d["Power"])
annotate_points(ax, data_2d["N"], data_2d["Speedup"], data_2d["Power"])
annotate_points(ax, data_3d["N"], data_3d["Speedup"], data_3d["Power"])

ax.axhline(1, linestyle="--", color="black")
ax.set_xlabel("Grid Size N")
ax.set_ylabel("Speedup (CPU / GPU)")
ax.set_title("GPU Speedup Across Dimensions")
ax.legend()

plt.tight_layout()
plt.savefig("speedup_all_dimensions.png", dpi=600)
plt.show()


# -------------------------------------------------------
# 3️⃣ Empirical Complexity (Slope Estimation)
# -------------------------------------------------------
def compute_slope(data, label):

    valid = data.dropna()

    x = np.log(valid["N"])
    y = np.log(valid["GPU"])

    slope = np.polyfit(x, y, 1)[0]
    print(f"Empirical GPU complexity for {label}: O(N^{slope:.2f})")


print("\n---- Empirical Complexity ----")
compute_slope(data_1d, "1D")
compute_slope(data_2d, "2D")
compute_slope(data_3d, "3D")


# -------------------------------------------------------
# 4️⃣ Efficiency Trend (Normalized Speedup)
# -------------------------------------------------------
fig, ax = plt.subplots(figsize=(10,7))

ax.plot(data_1d["Power"], data_1d["Speedup"], 'o-', label="1D")
ax.plot(data_2d["Power"], data_2d["Speedup"], 's-', label="2D")
ax.plot(data_3d["Power"], data_3d["Speedup"], '^-', label="3D")

annotate_points(ax, data_1d["Power"], data_1d["Speedup"], data_1d["Power"])
annotate_points(ax, data_2d["Power"], data_2d["Speedup"], data_2d["Power"])
annotate_points(ax, data_3d["Power"], data_3d["Speedup"], data_3d["Power"])

ax.set_xlabel("Power")
ax.set_ylabel("Speedup")
ax.set_title("Speedup vs Problem Growth")
ax.legend()

plt.tight_layout()
plt.savefig("efficiency_trend.png", dpi=600)
plt.show()



# print("\nBeautiful performance analysis completed 🚀")
