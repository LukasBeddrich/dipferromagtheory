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

nrhos = [61, 101, 141, 201]

dsfs = []
for nrho in nrhos:
    dsf = DynScalingFunc(**params)
    dsf.set("nrho", nrho)
    dsfs.append(dsf)

iterations = 50
save_the_arrays = np.zeros((len(nrhos), iterations, 2, len(params["qs"])))

t0 = time.time()
for i in range(iterations):
    for n, dsf in enumerate(dsfs):
        start = time.time()
#        save_the_arrays[n, i, 0], save_the_arrays[n, i, 1] = dsf.calc()
        finish = time.time()
        print(f"Iteration ({i+1}, {n+1}) took: {finish-start:.2f} seconds.")
print(f"Total calculation time: {finish-t0:.2f} seconds.")

#savedicts = {f'{nrho}' : save_the_arrays[ridx] for ridx, nrho in enumerate(nrhos)}
#np.savez(resdir + "/rho_points_study_results.npz", **savedicts)

storedarray = np.load(resdir + "/rho_points_study_results.npz")

plt.figure(figsize=(9,6))
for idx, nrho in enumerate(nrhos):
    plt.plot(dsfs[idx].x, storedarray[f"{nrho}"][-1, 0], ls="--", marker=".", label=f"n$_\\rho$ = {nrho}")
    #plt.plot(dsfs[idx].x, storedarray[f"{nrho}"][-1, 0] / storedarray[f"{nrho}"][-1, 0, -1], ls="--", marker=".", label=f"n$_\\rho$ = {nrho}")
    #plt.plot(dsfs[idx].x, dsfs[idx].gammaL / dsfs[idx].gammaL[-1], ls="--", marker=".", label=f"n$_\\rho$ = {nrho}")
plt.xlabel("$1 / q \\xi$", fontsize=16)
plt.ylabel("$\\gamma^L / \\gamma^L_0$", fontsize=16)
plt.legend()

plt.figure(figsize=(9,6))
for idx, nrho in enumerate(nrhos):
    plt.plot(dsfs[idx].x, storedarray[f"{nrho}"][-1, 1], ls="--", marker=".", label=f"n$_\\rho$ = {nrho}")
    #plt.plot(dsfs[idx].x, storedarray[f"{nrho}"][-1, 1] / storedarray[f"{nrho}"][-1, 1, -1], ls="--", marker=".", label=f"n$_\\rho$ = {nrho}")
    #plt.plot(dsfs[idx].x, dsfs[idx].gammaT / dsfs[idx].gammaT[-1], ls="--", marker=".", label=f"n$_\\rho$ = {nrho}")
plt.xlabel("$1 / q \\xi$", fontsize=16)
plt.ylabel("$\\gamma^T / \\gamma^T_0$", fontsize=16)
plt.legend()
plt.show()