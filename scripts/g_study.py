import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import time

from dipferromagtheory.linewidth.linewidth import DynScalingFunc
from dipferromagtheory import resdir

"""
Je größer g ist, umso mehr iterationen sind nötig?
"""

### Parameter set
gs = [0.0, 0.1, 1.0, 5.0, 10.0, 50.0, 100.0]

params = lambda g: dict(j=10.0, g=g, xi=1.0, qs=np.logspace(-1, 3, 161))

dsfs = []
for g in gs:
    dsf = DynScalingFunc(**params(g))
    dsfs.append(dsf)

iterations = 50
save_the_arrays = np.zeros((len(gs), iterations, 2, 161))

t0 = time.time()
for i in range(iterations):
    for n, dsf in enumerate(dsfs):
        start = time.time()
#        save_the_arrays[n, i, 0], save_the_arrays[n, i, 1] = dsf.calc()
        finish = time.time()
        print(f"Iteration ({i+1}, {n+1}) took: {finish-start:.2f} seconds.")
print(f"Total calculation time: {finish-t0:.2f} seconds.")

#savedicts = {f'{g}' : save_the_arrays[gidx] for gidx, g in enumerate(gs)}
#np.savez(resdir + "/g_study_results.npz", **savedicts)

storedarray = np.load(resdir + "/g_study_results.npz")

plt.figure(figsize=(9,6))
for idx, g in enumerate(gs):
    plt.plot(dsfs[idx].x, storedarray[f"{g}"][-1, 0]/storedarray[f"{g}"][-1, 0,-1], ls="--", marker=".", label=f"g = {g:.2f}")
#    plt.plot(dsfs[idx].x, dsfs[idx].gammaL / dsfs[idx].gammaL[-1], ls="--", marker=".", label=f"g = {g:.2f}")
#    plt.xlim((0,8))
#    plt.ylim((0.5, 3.0))
plt.xlabel("$1 / q \\xi$", fontsize=20)
plt.ylabel("$\\gamma^L / \\gamma^L_0$", fontsize=20)
plt.legend()

plt.figure(figsize=(9,6))
for idx, g in enumerate(gs):
    plt.plot(dsfs[idx].x, storedarray[f"{g}"][-1, 1]/storedarray[f"{g}"][-1, 1,-1], ls="--", marker=".", label=f"g = {g:.2f}")
#    plt.plot(dsfs[idx].x, dsfs[idx].gammaT / dsfs[idx].gammaT[-1], ls="--", marker=".", label=f"g = {g:.2f}")
#    plt.xlim((0,8))
#    plt.ylim((0.5, 3.0))
plt.xlabel("$1 / q \\xi$", fontsize=20)
plt.ylabel("$\\gamma^T / \\gamma^T_0$", fontsize=20)
plt.legend()
plt.show()