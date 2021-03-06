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
svs = [0.1, 1.0, 10.0, 100.0]

params = dict(
    j = 10.0,
    g = 5.0,
    xi = 1.0,
    qs = np.logspace(-1, 3, 161)
)

dsfs = []
for sv in svs:
    dsf = DynScalingFunc(**params)
    dsf.gammaL = np.ones(dsf.gammaL.shape) * sv
    dsf.gammaT = np.ones(dsf.gammaT.shape) * sv
    dsfs.append(dsf)

iterations = 50
save_the_arrays = np.zeros((len(svs), iterations, 2, 161))

t0 = time.time()
for i in range(iterations):
    for n, dsf in enumerate(dsfs):
        start = time.time()
#        save_the_arrays[n, i, 0], save_the_arrays[n, i, 1] = dsf.calc()
        finish = time.time()
        print(f"Iteration ({i+1}, {n+1}) took: {finish-start:.2f} seconds.")
print(f"Total calculation time: {finish-t0:.2f} seconds.")

#savedicts = {f'{sv}' : save_the_arrays[sidx] for sidx, sv in enumerate(svs)}
#np.savez(resdir + "/start_value_study_results.npz", **savedicts)

storedarray = np.load(resdir + "/start_value_study_results.npz")

plt.figure(figsize=(9,6))
for idx, sv in enumerate(svs):
#    plt.plot(dsfs[idx].x, dsfs[idx].gammaL / dsfs[idx].gammaL[-1], ls="--", marker=".", label=f"sv = {sv:.2f}")
    plt.plot(dsfs[idx].x, storedarray[f"{sv}"][-1, 0] / storedarray[f"{sv}"][-1, 0, -1], ls="--", marker=".", label=f"sv = {sv:.1f}")
plt.xlabel("$1 / q \\xi$", fontsize=20)
plt.ylabel("$\\gamma^L / \\gamma^L_0$", fontsize=20)
plt.legend()

plt.figure(figsize=(9,6))
for idx, sv in enumerate(svs):
#    plt.plot(dsfs[idx].x, dsfs[idx].gammaT / dsfs[idx].gammaT[-1], ls="--", marker=".", label=f"sv = {sv:.2f}")
    plt.plot(dsfs[idx].x, storedarray[f"{sv}"][-1, 1] / storedarray[f"{sv}"][-1, 1, -1], ls="--", marker=".", label=f"sv = {sv:.1f}")
plt.xlabel("$1 / q \\xi$", fontsize=20)
plt.ylabel("$\\gamma^T / \\gamma^T_0$", fontsize=20)
plt.legend()
plt.show()