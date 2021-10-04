import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import time

from dipferromagtheory.linewidth.linewidth import DynScalingFunc
from dipferromagtheory import resdir

"""

"""

### Parameter set

params = dict(
    j = 10.0,
    g = 5.0,
    xi = 1.0,
    qs = np.logspace(-1, 3, 161)
)

dsf = DynScalingFunc(**params)

iterations = [10, 20, 50, 100, 200]
save_the_arrays = np.zeros((max(iterations), 2, 161))

t0 = time.time()
for i in range(max(iterations)):
    start = time.time()
    save_the_arrays[i, 0], save_the_arrays[i, 1] = dsf.calc()
    finish = time.time()
    print(f"Iteration ({i+1}) took: {finish-start:.2f} seconds.")
print(f"Total calculation time: {finish-t0:.2f} seconds.")

savedict = {f'{max(iterations)}' : save_the_arrays}
np.savez(resdir + "/iterations_study_results.npz", **savedict)

storedarray = np.load(resdir + "/iterations_study_results.npz")

plt.figure(figsize=(9,6))
for n in range(max(iterations)):
    if n+1 in iterations:
        plt.plot(dsf.x, storedarray[f"{max(iterations)}"][n, 0], ls="--", marker=".", label=f"n = {n+1}")
#       plt.plot(dsf.x, dsfs[idx].gammaL, ls="--", marker=".", label=f"n = {n:.2f}")
plt.legend()

plt.figure(figsize=(9,6))
for n in range(max(iterations)):
    if n+1 in iterations:
        plt.plot(dsf.x, storedarray[f"{max(iterations)}"][n, 1], ls="--", marker=".", label=f"n = {n+1}")
#    plt.plot(dsfs[idx].x, dsfs[idx].gammaT, ls="--", marker=".", label=f"n = {n:.2f}")
plt.legend()
plt.show()