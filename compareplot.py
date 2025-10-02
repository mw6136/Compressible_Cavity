import numpy as np
import glob
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable

field_key = {
    "Density": "Density",
    "VelX1": "X-Velocity",
    "VelX2": "Y-Velocity",
    "Press": "Pressure",
    "Temp": "Temperature",
    "VelMag": "Velocity Magnitude",
    "dpdx": "dpdx",
    "dpdy": "dpdy",
}

dt = 0.1

Time = 9

ds_08 = h5py.File(f"comparison_0.8/test.output.{str(Time).zfill(4)}.h5", "r")
ds_05 = h5py.File(f"comparison_0.5/test.output.{str(Time).zfill(4)}.h5", "r")
ds_01 = h5py.File(f"comparison_0.1/test.output.{str(Time).zfill(4)}.h5", "r")
ds_25 = h5py.File(f"comparison_0.025/test.output.{str(Time).zfill(4)}.h5", "r")

gmax, gmin = 0.0, 1e20

N = 64
field = "VelMag"

for i, ds in enumerate([ds_25, ds_01, ds_05, ds_08]):

    u = np.rot90(ds["VelX1"][:].reshape(N + 2, N + 2)[1:-1, 1:-1], 3)
    v = np.rot90(ds["VelX2"][:].reshape(N + 2, N + 2)[1:-1, 1:-1], 3)

    plotarray = np.sqrt(u**2 + v**2) / np.sqrt(u**2 + v**2).max()

    pmax = np.amax(plotarray)
    pmin = np.amin(plotarray)
    if pmax > gmax:
        gmax = pmax
    else:
        gmax = gmax

    if pmin < gmin:
        gmin = pmin
    else:
        gmin = gmin


fig, ax = plt.subplots(1, 4, figsize=(18, 6), sharey=True)

Mas = [0.025, 0.1, 0.5, 0.8]

for i, ds in enumerate([ds_25, ds_01, ds_05, ds_08]):

    u = np.rot90(ds["VelX1"][:].reshape(N + 2, N + 2)[1:-1, 1:-1], 3)
    v = np.rot90(ds["VelX2"][:].reshape(N + 2, N + 2)[1:-1, 1:-1], 3)

    plotarray = np.sqrt(u**2 + v**2) / np.sqrt(u**2 + v**2).max()

    x, y = np.linspace(0, 1, len(plotarray)), np.linspace(0, 1, len(plotarray))
    X, Y = np.meshgrid(x, y)

    im = ax[i].pcolormesh(X, Y, plotarray, cmap="jet", vmax=gmax, vmin=gmin)
    ax[i].set_title(f"$\mathrm{{Ma}} = ${Mas[i]}")
    if i == 3:
        divider = make_axes_locatable(ax[i])
        cax = divider.append_axes("right", size="5%", pad=0.02)
        fig.add_axes(cax)
        fig.colorbar(im, cax=cax)

fig.suptitle(
    rf"Normalized {field_key[field]}, $\omega t = ${round(float(Time) * float(dt),7)}", fontsize=16
)
fig.tight_layout()

plt.savefig(f"Ma_comparison.png", bbox_inches="tight", dpi=300)
plt.show()
plt.close()
