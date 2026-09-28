"""
Identify the terminase large subunit (TerL) protein in the Bakta annotation
by its characteristic Pfam domain signature(s), and write it out as a
single-sequence FASTA for use as the tblastn query in the phylogenetics rule.

Injected by Snakemake: snakemake.input.done, snakemake.output.faa,
snakemake.params.annotation_dir, snakemake.params.pfam_domains
"""

import json
from pathlib import Path

pfam_domains = set(snakemake.params.pfam_domains)
out_faa = snakemake.output.faa
bakta_json = Path(snakemake.input.json)

with open(bakta_json) as f:
    data = json.load(f)

hit_feature = None
for feature in data.get("features", []):
    # Collect all string values that might carry a Pfam accession, from
    # any of Bakta's cross-reference / domain-search fields.
    blob = json.dumps(feature)
    if any(domain in blob for domain in pfam_domains):
        hit_feature = feature
        break

if hit_feature is None:
    raise RuntimeError(
        f"No feature matched Pfam domains {pfam_domains} in {bakta_json}. "
        "Inspect the Bakta JSON manually and adjust terl_pfam_domains in "
        "config.yaml, or identify TerL by its expected genomic position "
        "(immediately downstream of the terminase small subunit gene)."
    )

locus_tag = hit_feature.get("locus", hit_feature.get("locus_tag", "TerL"))
aa_seq = hit_feature.get("aa") or hit_feature.get("sequence")
if not aa_seq:
    raise RuntimeError(f"Matched feature {locus_tag} has no translated sequence field.")

with open(out_faa, "w") as out:
    out.write(f">{locus_tag} terminase_large_subunit\n{aa_seq}\n")

print(f"Identified TerL as {locus_tag} ({len(aa_seq)} aa) -> {out_faa}")
