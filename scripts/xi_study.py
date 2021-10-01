import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import time

from dipferromagtheory.linewidth.linewidth import DynScalingFunc
from dipferromagtheory import resdir

"""

"""

### Parameter set
xis = [0.1, 1.0, 10.0]

params = lambda xi: dict(j=10.0, g=10.0, xi=xi, qs=np.logspace(-1, 3, 101)/xi)

dsfs = []
for xi in xis:
    dsf = DynScalingFunc(**params(xi))
    dsfs.append(dsf)

iterations = 20
save_the_arrays = np.zeros((len(xis), iterations, 2, 101))

t0 = time.time()
for i in range(iterations):
    for n, dsf in enumerate(dsfs):
        start = time.time()
        save_the_arrays[n, i, 0], save_the_arrays[n, i, 1] = dsf.calc()
        finish = time.time()
        print(f"Iteration ({i+1}, {n+1}) took: {finish-start:.2f} seconds.")
print(f"Total calculation time: {finish-t0:.2f} seconds.")

savedicts = {f'{xi}' : save_the_arrays[xidx] for xidx, xi in enumerate(xis)}
np.savez(resdir + "/xi_study_results.npz", **savedicts)

storedarray = np.load(resdir + "/xi_study_results.npz")

plt.figure(figsize=(9,6))
for idx, xi in enumerate(xis):
#    plt.plot(dsfs[idx].x, storedarray[f"{xi}"][-1, 0], ls="--", marker=".", label=f"$\\xi$ = {xi:.2f}")
    plt.plot(dsfs[idx].x, dsfs[idx].gammaL, ls="--", marker=".", label=f"$\\xi$ = {xi:.2f}")
plt.legend()

plt.figure(figsize=(9,6))
for idx, xi in enumerate(xis):
#    plt.plot(dsfs[idx].x, storedarray[f"{xi}"][-1, 1], ls="--", marker=".", label=f"$\\xi$ = {xi:.2f}")
    plt.plot(dsfs[idx].x, dsfs[idx].gammaT, ls="--", marker=".", label=f"$\\xi$ = {xi:.2f}")
plt.legend()
plt.show()