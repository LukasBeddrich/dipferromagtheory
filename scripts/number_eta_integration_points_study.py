import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import time

from dipferromagtheory.linewidth.linewidth import DynScalingFunc
from dipferromagtheory import resdir

### Parameter set
params = dict(
    j = 10.0,
    g = 5.0,
    xi = 1.0,
    qs = np.logspace(-1, 3, 161)
)

netas = [11, 21, 41, 81]

dsfs = []
for neta in netas:
    dsf = DynScalingFunc(**params)
    dsf.set("neta", neta)
    dsfs.append(dsf)

iterations = 50
save_the_arrays = np.zeros((len(netas), iterations, 2, len(params["qs"])))

t0 = time.time()
for i in range(iterations):
    for n, dsf in enumerate(dsfs):
        start = time.time()
        save_the_arrays[n, i, 0], save_the_arrays[n, i, 1] = dsf.calc()
        finish = time.time()
        print(f"Iteration ({i+1}, {n+1}) took: {finish-start:.2f} seconds.")
print(f"Total calculation time: {finish-t0:.2f} seconds.")

savedicts = {f'{neta}' : save_the_arrays[ridx] for ridx, neta in enumerate(netas)}
np.savez(resdir + "/eta_points_study_results.npz", **savedicts)

#storedarray = np.load(resdir + "/eta_points_study_results.npz")

plt.figure(figsize=(9,6))
for idx, neta in enumerate(netas):
    #plt.plot(dsfs[idx].x, storedarray[f"{rul}_{neta}"][-1, 0] / storedarray[f"{rul}_{neta}"][-1, 0, -1], ls="--", marker=".", label=f"n$_\\eta$ = {neta}")
    plt.plot(dsfs[idx].x, dsfs[idx].gammaL / dsfs[idx].gammaL[-1], ls="--", marker=".", label=f"n$_\\eta$ = {neta}")
plt.xlabel("$1 / q \\xi$", fontsize=16)
plt.ylabel("$\\gamma^L / \\gamma^L_0$", fontsize=16)
plt.legend()

plt.figure(figsize=(9,6))
for idx, neta in enumerate(netas):
    #plt.plot(dsfs[idx].x, storedarray[f"{rul}_{neta}"][-1, 1] / storedarray[f"{rul}_{neta}"][-1, 0, -1], ls="--", marker=".", label=f"n$_\\eta$ = {neta}")
    plt.plot(dsfs[idx].x, dsfs[idx].gammaT / dsfs[idx].gammaT[-1], ls="--", marker=".", label=f"n$_\\eta$ = {neta}")
plt.xlabel("$1 / q \\xi$", fontsize=16)
plt.ylabel("$\\gamma^T / \\gamma^T_0$", fontsize=16)
plt.legend()
plt.show()