"""
Summarize per-base depth into fixed-size windows and flag windows whose mean
depth deviates sharply from the genome-wide mean -- a signature of either a
chimeric assembly (two different genomes stitched together) or a
collapsed/duplicated repeat region.

Injected by Snakemake: snakemake.input.raw, snakemake.output.summary,
snakemake.output.flagged, snakemake.params.window
"""

import numpy as np

raw_path = snakemake.input.raw
summary_path = snakemake.output.summary
flagged_path = snakemake.output.flagged
window = snakemake.params.window

positions, depths = [], []
with open(raw_path) as f:
    for line in f:
        _, pos, d = line.split()
        positions.append(int(pos))
        depths.append(int(d))

depths = np.array(depths)
n_windows = len(depths) // window

means = []
with open(summary_path, "w") as out:
    out.write("window_start\tmean_depth\n")
    for i in range(n_windows):
        seg = depths[i * window:(i + 1) * window]
        m = seg.mean()
        means.append(m)
        out.write(f"{i * window}\t{m:.1f}\n")

means = np.array(means)
genome_mean = means.mean()
genome_std = means.std()

with open(flagged_path, "w") as out:
    out.write("window_start\tmean_depth\tgenome_mean\tgenome_std\treason\n")
    for i, m in enumerate(means):
        if m < 0.2 * genome_mean:
            out.write(f"{i*window}\t{m:.1f}\t{genome_mean:.1f}\t{genome_std:.1f}\tlow_depth_possible_chimera_junction\n")
        elif m > 3 * genome_mean:
            out.write(f"{i*window}\t{m:.1f}\t{genome_mean:.1f}\t{genome_std:.1f}\thigh_depth_possible_collapsed_repeat\n")

print(f"Genome-wide mean depth: {genome_mean:.1f} (SD {genome_std:.1f})")
print(f"Windows evaluated: {n_windows}")
n_flagged = sum(1 for m in means if m < 0.2 * genome_mean or m > 3 * genome_mean)
print(f"Flagged windows: {n_flagged}")
if n_flagged == 0:
    print("No abrupt depth discontinuities detected -- no evidence of "
          "chimeric assembly or collapsed repeats from this check.")
else:
    print("WARNING: flagged windows found -- inspect flagged_windows.tsv "
          "and consider whether the assembly may be chimeric or contain "
          "a collapsed repeat.")
