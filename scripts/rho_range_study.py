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

# try to keep number of points for each order of magnitude constant
rho_upper_limits = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
nrhos = [61, 71, 81, 91, 101, 111, 121, 131, 141]

dsfs = []
for rul, nrho in zip(rho_upper_limits, nrhos):
    dsf = DynScalingFunc(**params)
    dsf.set("rho_upper_limit", rul)
    dsf.set("nrho", nrho)
    dsfs.append(dsf)

iterations = 50
save_the_arrays = np.zeros((len(rho_upper_limits), iterations, 2, len(params["qs"])))

t0 = time.time()
for i in range(iterations):
    for n, dsf in enumerate(dsfs):
        start = time.time()
#        save_the_arrays[n, i, 0], save_the_arrays[n, i, 1] = dsf.calc()
        finish = time.time()
        print(f"Iteration ({i+1}, {n+1}) took: {finish-start:.2f} seconds.")
print(f"Total calculation time: {finish-t0:.2f} seconds.")

#savedicts = {f'{rul}_{nrho}' : save_the_arrays[ridx] for ridx, (rul, nrho) in enumerate(zip(rho_upper_limits, nrhos))}
#np.savez(resdir + "/rho_range_study_results.npz", **savedicts)

storedarray = np.load(resdir + "/rho_range_study_results.npz")

plt.figure(figsize=(9,6))
for idx, (rul, nrho) in enumerate(zip(rho_upper_limits, nrhos)):
    plt.plot(dsfs[idx].x, storedarray[f"{rul}_{nrho}"][-1, 0] / storedarray[f"{rul}_{nrho}"][-1, 0, -1], ls="--", marker=".", label=f"$\\rho_{'{max}'}$ = {rul:.1f} | n$_\\rho$ = {nrho}")
    #plt.plot(dsfs[idx].x, dsfs[idx].gammaL, ls="--", marker=".", label=f"$\\rho_{'{max}'}$ = {rul:.2f} | n$_\\rho$ = {nrho}")
plt.xlabel("$1 / q \xi$", fontsize=20)
plt.ylabel("$\gamma^L / \gamma^L_0$", fontsize=20)
plt.legend()

plt.figure(figsize=(9,6))
for idx, (rul, nrho) in enumerate(zip(rho_upper_limits, nrhos)):
    plt.plot(dsfs[idx].x, storedarray[f"{rul}_{nrho}"][-1, 1] / storedarray[f"{rul}_{nrho}"][-1, 0, -1], ls="--", marker=".", label=f"$\\rho_{'{max}'}$ = {rul:.1f} | n$_\\rho$ = {nrho}")
    #plt.plot(dsfs[idx].x, dsfs[idx].gammaT, ls="--", marker=".", label=f"$\\rho_{'{max}'}$ = {rul:.2f} | n$_\\rho$ = {nrho}")
plt.xlabel("$1 / q \xi$", fontsize=20)
plt.ylabel("$\gamma^T / \gamma^T_0$", fontsize=20)
plt.legend()
plt.show()