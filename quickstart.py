"""RSE vs sketch size (m) for different sketch families."""
import matplotlib.pyplot as plt
from weighted_cardinality_estimation import FastExpSketch, QSketch, kQSketchShifted, HyperLogLog, MemoryFlag, stat

ms, n_reps, N = list(range(30, 61)), 10_000, 500
elements, weights = stat.weighted_stream(n=N, total_weight=N, dist="zipf", seed=0)

sketches = {
    "FastExpSketch (64-bit)": lambda m, s: FastExpSketch(m, seed=s),
    "QSketch (8-bit)": lambda m, s: QSketch(m, seed=s, amount_bits=8),
    "kQSketchShifted (5-bit, k=2)": lambda m, s: kQSketchShifted(m, seed=s, amount_bits=5, logarithm_base=2.0),
    "HyperLogLog": lambda m, s: HyperLogLog(m, seed=s),
}

results = {name: [] for name in sketches}
memory = {name: [] for name in sketches}

for m in ms:
    for name, make in sketches.items():
        errs = []
        for rep in range(n_reps):
            sk = make(m, rep)
            if "HyperLogLog" in name:
                sk.add_many(elements)
            else:
                sk.add_many(elements, weights)
            errs.append(sk.estimate() / N - 1)
        results[name].append(stat.compute_rse(errs))
        memory[name].append(make(m, 0).memory_usage(MemoryFlag.TOTAL))

import numpy as np

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

for name, rses in results.items():
    ax1.plot(ms, rses, marker="o", markersize=3, label=name)
ax1.plot(ms, 1 / np.sqrt(np.array(ms) - 2), "k--", label=r"$1/\sqrt{m-2}$")
ax1.set_xlabel("Sketch size (m)")
ax1.set_ylabel("Relative Standard Error")
ax1.legend()
ax1.grid(True, alpha=0.3)

for name, mem in memory.items():
    ax2.plot(mem, results[name], marker="o", markersize=3, label=name)
ax2.set_xlabel("Memory (bytes)")
ax2.set_ylabel("Relative Standard Error")
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("quickstart_plot.png", dpi=150)
