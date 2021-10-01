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
js = [0.0, 0.1, 1.0, 10.0, 100.0]

params = lambda j: dict(j=j, g=10.0, xi=1.0, qs=np.logspace(-1, 3, 101))

dsfs = []
for j in js:
    dsf = DynScalingFunc(**params(j))
    dsfs.append(dsf)

iterations = 20
save_the_arrays = np.zeros((len(js), iterations, 2, 101))

t0 = time.time()
for i in range(iterations):
    for n, dsf in enumerate(dsfs):
        start = time.time()
        save_the_arrays[n, i, 0], save_the_arrays[n, i, 1] = dsf.calc()
        finish = time.time()
        print(f"Iteration ({i+1}, {n+1}) took: {finish-start:.2f} seconds.")
print(f"Total calculation time: {finish-t0:.2f} seconds.")

savedicts = {f'{j}' : save_the_arrays[jidx] for jidx, j in enumerate(js)}
np.savez(resdir + "/j_study_results.npz", **savedicts)

storedarray = np.load(resdir + "/j_study_results.npz")

plt.figure(figsize=(9,6))
for idx, j in enumerate(js):
#    plt.plot(dsfs[idx].x, storedarray[f"{j}"][-1, 0], ls="--", marker=".", label=f"j = {j:.2f}")
    plt.plot(dsfs[idx].x, dsfs[idx].gammaL, ls="--", marker=".", label=f"j = {j:.2f}")
plt.legend()

plt.figure(figsize=(9,6))
for idx, j in enumerate(js):
#    plt.plot(dsfs[idx].x, storedarray[f"{j}"][-1, 1], ls="--", marker=".", label=f"j = {j:.2f}")
    plt.plot(dsfs[idx].x, dsfs[idx].gammaT, ls="--", marker=".", label=f"j = {j:.2f}")
plt.legend()
plt.show()