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

rho_upper_limits = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]

dsfs = []
for rul in rho_upper_limits:
    dsf = DynScalingFunc(**params)
    dsf.set("rho_upper_limit", rul)
    dsfs.append(dsf)

iterations = 10
save_the_arrays = np.zeros((len(rho_upper_limits), iterations, 2, len(params["qs"])))

t0 = time.time()
for i in range(iterations):
    for n, dsf in enumerate(dsfs):
        start = time.time()
#        save_the_arrays[n, i, 0], save_the_arrays[n, i, 1] = dsf.calc()
        finish = time.time()
        print(f"Iteration ({i+1}, {n+1}) took: {finish-start:.2f} seconds.")
print(f"Total calculation time: {finish-t0:.2f} seconds.")

#savedicts = {f'{rul}' : save_the_arrays[ridx] for ridx, rul in enumerate(rho_upper_limits)}
#np.savez(resdir + "/rho_range_study_results.npz", **savedicts)

storedarray = np.load(resdir + "/rho_range_study_results.npz")

plt.figure(figsize=(9,6))
for idx, rul in enumerate(rho_upper_limits):
    plt.plot(dsfs[idx].x, storedarray[f"{rul}"][-1, 0], ls="--", marker=".", label=f"$\\rho_{'{max}'}$ = {rul:.2f}")
    #plt.plot(dsfs[idx].x, dsfs[idx].gammaL, ls="--", marker=".", label=f"$\\rho_m$ = {np.exp(rul):.2f}")
plt.legend()

plt.figure(figsize=(9,6))
for idx, rul in enumerate(rho_upper_limits):
    plt.plot(dsfs[idx].x, storedarray[f"{rul}"][-1, 1], ls="--", marker=".", label=f"$\\rho_m$ = {rul:.2f}")
    #plt.plot(dsfs[idx].x, dsfs[idx].gammaT, ls="--", marker=".", label=f"$\\rho_m$ = {np.exp(rul):.2f}")
plt.legend()
plt.show()