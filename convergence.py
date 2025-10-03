import numpy as np
import matplotlib.pyplot as plt
import glob
from tqdm import tqdm
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse

###################################### Spatial Convergence #################################################
Time = 42

ds_16 = h5py.File(f"data/convergence_16x16/test.output.{str(Time).zfill(4)}.h5", "r")
ds_16 = np.rot90(ds_16["VelX1"][:].reshape(16 + 2, 16 + 2)[1:-1, 1:-1], 3)

ds_32 = h5py.File(f"data/convergence_32x32/test.output.{str(Time).zfill(4)}.h5", "r")
ds_32 = np.rot90(ds_32["VelX1"][:].reshape(32 + 2, 32 + 2)[1:-1, 1:-1], 3)

ds_64 = h5py.File(f"data/convergence_64x64/test.output.{str(Time).zfill(4)}.h5", "r")
ds_64 = np.rot90(ds_64["VelX1"][:].reshape(64 + 2, 64 + 2)[1:-1, 1:-1], 3)

ds_128 = h5py.File(f"data/convergence_128x128/test.output.{str(Time).zfill(4)}.h5", "r")
ds_128 = np.rot90(ds_128["VelX1"][:].reshape(128 + 2, 128 + 2)[1:-1, 1:-1], 3)

eps_128_64 = np.sqrt((ds_128[128 // 2, 128 // 2] - ds_64[64 // 2, 64 // 2]) ** 2)
eps_128_32 = np.sqrt((ds_128[128 // 2, 128 // 2] - ds_32[32 // 2, 32 // 2]) ** 2)
eps_128_16 = np.sqrt((ds_128[128 // 2, 128 // 2] - ds_16[16 // 2, 16 // 2]) ** 2)

dxs = [1 / 128, 1 / 64, 1 / 32]
eps = [eps_128_64, eps_128_32, eps_128_16]

plt.style.use("classic")
fig = plt.figure(figsize=(8, 8))

plt.plot(dxs, eps, "--", marker="o", markersize=12, linewidth=3, color="blue")
plt.plot(
    0.2 * (np.logspace(-1.8, -1.1, 100)),
    10 * (np.logspace(-1.8, -1.1, 100)) ** (2),
    "-",
    linewidth=3,
    color="black",
    label="$\propto (\Delta x)^{2}$",
)

plt.yscale("log")
plt.xscale("log")

plt.xlabel(r"Grid Spacing $\Delta x$", fontsize=16)
plt.ylabel(
    r"Norm Error $||\epsilon||_{L^2} = \sqrt{\frac{1}{N^2} \sum_{1 \leq i,j \leq N} [u \, (x_i, y_j) - u^*(x_i, y_j)]^2}$",
#    r"Norm Error $||\epsilon||_{L^2} = \sqrt{\frac{1}{N^2} \sum_{1 \leq j \leq N} [u \, (x_j, t) - u^*(x_j, t)]^2}$",
    fontsize=16,
)
plt.title(r"Spatial Convergence $\mathrm{Ma} = 0.025$", fontsize=16)

plt.minorticks_on()
plt.legend(loc="upper left")
plt.grid()

plt.savefig("spatial_convergence.png", bbox_inches="tight", dpi=300)
plt.show()
plt.close()

#print("slope = ",(eps[-1] - eps[0])/(dxs[-1] - dxs[0]))

###################################### Temporal Convergence ############################################

Time = 41

ds_01 = h5py.File(f"data/time_conv_01/test.output.{str(Time).zfill(4)}.h5", "r")
ds_01 = np.rot90(ds_01["VelX1"][:].reshape(32 + 2, 32 + 2)[1:-1, 1:-1], 3)

ds_005 = h5py.File(f"data/time_conv_005/test.output.{str(Time).zfill(4)}.h5", "r")
ds_005 = np.rot90(ds_005["VelX1"][:].reshape(32 + 2, 32 + 2)[1:-1, 1:-1], 3)

ds_001 = h5py.File(f"data/time_conv_001/test.output.{str(Time).zfill(4)}.h5", "r")
ds_001 = np.rot90(ds_001["VelX1"][:].reshape(32 + 2, 32 + 2)[1:-1, 1:-1], 3)

ds_0005 = h5py.File(f"data/time_conv_0005/test.output.{str(Time).zfill(4)}.h5", "r")
ds_0005 = np.rot90(ds_0005["VelX1"][:].reshape(32 + 2, 32 + 2)[1:-1, 1:-1], 3)

eps_0005_001 = np.sqrt((ds_0005[32 // 2, 32 // 2] - ds_001[32 // 2, 32 // 2]) ** 2)
eps_0005_005 = np.sqrt((ds_0005[32 // 2, 32 // 2] - ds_005[32 // 2, 32 // 2]) ** 2)
eps_0005_01 = np.sqrt((ds_0005[32 // 2, 32 // 2] - ds_01[32 // 2, 32 // 2]) ** 2)

dts = [4.39067e-8, 8.78134e-8, 4.39067e-7]
eps = [eps_0005_001, eps_0005_005, eps_0005_01]

plt.style.use("classic")
fig = plt.figure(figsize=(8, 8))

plt.plot(dts, eps, "--", marker="o", markersize=12, linewidth=3, color="red")
plt.plot(
    np.logspace(-7.5, -6.5, 100),
    20 * (np.logspace(-3, -1.5, 100)),
    "-",
    linewidth=3,
    color="black",
    label="$\propto \Delta t$",
)
plt.plot(
    np.logspace(-7.5, -6.5, 100),
    10 * (np.logspace(-7.5, -6.5, 100)) ** (1 / 2),
    "--",
    linewidth=3,
    color="black",
    label="$\propto (\Delta t)^{1/2}$",
)

plt.yscale("log")
plt.xscale("log")

plt.xlabel(r"Timestep $\Delta t$", fontsize=16)
plt.ylabel(
    r"Norm Error $||\epsilon||_{L^2} = \sqrt{\frac{1}{n^2} \sum_{1 \leq t_n \leq n} [u \, (t_n) - u^*(t_n)]^2}$",
    fontsize=16,
)
plt.title(r"Temporal Convergence $\mathrm{Ma} = 0.025$", fontsize=16)

plt.minorticks_on()
plt.legend(loc="upper left")
plt.grid()

plt.savefig("temporal_convergence.png", bbox_inches="tight", dpi=300)
plt.show()
plt.close()

#print("slope = ",(eps[-1] - eps[0])/(dts[-1] - dts[0]))
