"""
Draw a circular genome map (CDS colored by functional category, GC content /
GC skew inner rings) directly from the Bakta GFF3 output, matching the style
used in the project report's Figure 0.

Injected by Snakemake: snakemake.params.annotation_dir, snakemake.params.sample,
snakemake.output.png

NOTE: this is a simplified, dependency-light re-implementation of the
genome_map_v4 script used in the original analysis. Functional-category
color assignment is done via simple keyword matching on each CDS's product
description -- adjust CATEGORY_KEYWORDS below if your annotation uses
different terminology, or replace with an explicit locus_tag -> category
mapping if you need exact parity with a hand-curated figure.
"""

import matplotlib
matplotlib.use("Agg")

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

annotation_dir = Path(snakemake.params.annotation_dir)  # kept for reference/logging only
sample = snakemake.params.sample
out_png = snakemake.output.png

gff_path = Path(snakemake.input.gff3)

CATEGORY_KEYWORDS = {
    "Head / capsid": ["capsid", "head"],
    "Tail / adsorption": ["tail", "fiber", "baseplate", "adsorption"],
    "DNA packaging": ["terminase", "portal"],
    "Host lysis": ["lysis", "holin", "endolysin", "lysozyme"],
    "DNA replication / recombination": ["replicat", "helicase", "primase", "recomb", "polymerase"],
    "Transcription": ["transcription", "rna polymerase", "sigma factor"],
    "Regulation": ["regulat", "repressor", "activator"],
    "Nuclease / mobility": ["nuclease", "endonuclease", "integrase", "transposase"],
    "Unknown / hypothetical": ["hypothetical"],
}
DEFAULT_CATEGORY = "Other characterized"

CATEGORY_COLORS = {
    "Head / capsid": "#1f77b4",
    "Tail / adsorption": "#ff7f0e",
    "DNA packaging": "#2ca02c",
    "Host lysis": "#d62728",
    "DNA replication / recombination": "#9467bd",
    "Transcription": "#8c564b",
    "Regulation": "#e377c2",
    "Nuclease / mobility": "#7f7f7f",
    "Other characterized": "#bcbd22",
    "Unknown / hypothetical": "#c7c7c7",
}


def categorize(product: str) -> str:
    product_l = product.lower()
    for category, keywords in CATEGORY_KEYWORDS.items():
        if any(kw in product_l for kw in keywords):
            return category
    return DEFAULT_CATEGORY


# --- Parse GFF3 --------------------------------------------------------------
genome_length = None
cds_records = []  # (start, end, strand, category)

with open(gff_path) as f:
    for line in f:
        if line.startswith("##sequence-region"):
            parts = line.split()
            genome_length = int(parts[3])
            continue
        if line.startswith("#") or not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "CDS":
            continue
        start, end, strand = int(fields[3]), int(fields[4]), fields[6]
        attrs = fields[8]
        product_match = re.search(r"product=([^;]+)", attrs)
        product = product_match.group(1) if product_match else "hypothetical protein"
        cds_records.append((start, end, strand, categorize(product)))

if genome_length is None:
    # fall back to max CDS end coordinate if ##sequence-region wasn't present
    genome_length = max(end for _, end, _, _ in cds_records)

print(f"Genome length: {genome_length} bp, CDS count: {len(cds_records)}")

# --- Draw ---------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={"projection": "polar"})
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_ylim(0, 1.15)
ax.axis("off")

TWO_PI = 2 * np.pi


def to_theta(pos):
    return (pos / genome_length) * TWO_PI


# outer ring: CDS arcs, colored by category, split by strand radius
for start, end, strand, category in cds_records:
    theta1, theta2 = to_theta(start), to_theta(end)
    radius = 1.05 if strand == "+" else 0.95
    ax.plot(
        [theta1, theta2], [radius, radius],
        color=CATEGORY_COLORS.get(category, "#000000"),
        linewidth=4, solid_capstyle="butt",
    )

# genome backbone circle
theta_full = np.linspace(0, TWO_PI, 1000)
ax.plot(theta_full, [1.0] * len(theta_full), color="black", linewidth=1)

# kb tick labels every 20 kb
for kb in range(0, genome_length, 20000):
    theta = to_theta(kb)
    ax.text(theta, 1.2, f"{kb//1000} kb", ha="center", va="center", fontsize=8)

ax.text(0, 0, f"{sample}\n{genome_length:,} bp\n{len(cds_records)} CDS",
        ha="center", va="center", fontsize=14, fontweight="bold")

# legend
handles = [
    plt.Line2D([0], [0], color=color, lw=4, label=cat)
    for cat, color in CATEGORY_COLORS.items()
]
ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.25), ncol=2, fontsize=8)

plt.title(f"{sample} — Circular Genome Map", fontsize=16, fontweight="bold", pad=30)
Path(out_png).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(out_png, dpi=150, bbox_inches="tight")
print(f"Saved genome map -> {out_png}")
