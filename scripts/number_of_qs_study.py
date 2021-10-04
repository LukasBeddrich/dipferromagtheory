import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import time

from dipferromagtheory.linewidth.linewidth import DynScalingFunc
from dipferromagtheory import resdir

"""

"""

### Parameter set
nqs = [101, 151, 201, 401]

params = lambda n: dict(j=10.0, g=10.0, xi=1.0, qs=np.logspace(-1, 3, n))

dsfs = []
for n in nqs:
    dsf = DynScalingFunc(**params(n))
    dsfs.append(dsf)

iterations = 20
save_the_arrays = np.zeros((len(nqs), iterations, 2, max(nqs)))

t0 = time.time()
for i in range(iterations):
    for j, dsf in enumerate(dsfs):
        start = time.time()
#        save_the_arrays[j, i, 0, :nqs[j]], save_the_arrays[j, i, 1, :nqs[j]] = dsf.calc()
        finish = time.time()
#        print(f"Iteration ({i+1}, {j+1}) took: {finish-start:.2f} seconds.")
#print(f"Total calculation time: {finish-t0:.2f} seconds.")

#savedicts = {f'{n}' : save_the_arrays[nidx] for nidx, n in enumerate(nqs)}
#np.savez(resdir + "/number_of_qs_study_results.npz", **savedicts)

storedarray = np.load(resdir + "/number_of_qs_study_results.npz")

plt.figure(figsize=(9,6))
for idx, n in enumerate(nqs):
    plt.plot(dsfs[idx].x, storedarray[f"{n}"][-1, 0, :n], ls="--", marker=".", label=f"n = {n}")
#    plt.plot(dsfs[idx].x, dsfs[idx].gammaL, ls="--", marker=".", label=f"n = {n:.2f}")
plt.legend()

plt.figure(figsize=(9,6))
for idx, n in enumerate(nqs):
    plt.plot(dsfs[idx].x, storedarray[f"{n}"][-1, 1,:n], ls="--", marker=".", label=f"n = {n}")
#    plt.plot(dsfs[idx].x, dsfs[idx].gammaT, ls="--", marker=".", label=f"n = {n:.2f}")
plt.legend()
plt.show()