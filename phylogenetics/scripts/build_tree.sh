#!/bin/bash
# Terminase large subunit phylogeny across Krischvirus genomes
# Usage: run from a directory containing:
#   - a query fasta (terminase_large_query.fasta) with the seed protein sequence
#   - a phylo_genomes/ folder with one .fasta genome per taxon
set -e

QUERY="terminase_large_query.fasta"
GENOME_DIR="phylo_genomes"
ANALYSIS_DIR="phylo_analysis"
mkdir -p "$ANALYSIS_DIR"

echo "== tblastn: locating terminase large subunit homolog in each genome =="
for f in "$GENOME_DIR"/*.fasta; do
    name=$(basename "$f" .fasta)
    makeblastdb -in "$f" -dbtype nucl -out "$ANALYSIS_DIR/db_${name}" > /dev/null 2>&1
    tblastn -query "$QUERY" -db "$ANALYSIS_DIR/db_${name}" \
        -outfmt "6 qseqid sseqid pident length evalue bitscore sstart send sframe" \
        -evalue 1e-10 -max_target_seqs 1 > "$ANALYSIS_DIR/${name}_terminase_hit.tsv"
done

echo "== Extracting full nucleotide span + translating =="
python3 - <<PYEOF
import os
from Bio import SeqIO
from Bio.Seq import Seq

genome_dir = "$GENOME_DIR"
hit_dir = "$ANALYSIS_DIR"
results_nt, results_aa = [], []

for fname in os.listdir(hit_dir):
    if not fname.endswith("_terminase_hit.tsv"):
        continue
    name = fname.replace("_terminase_hit.tsv", "")
    coords, strand = [], None
    with open(os.path.join(hit_dir, fname)) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 9: continue
            sstart, send = int(parts[6]), int(parts[7])
            frame = int(parts[8])
            coords.append((sstart, send))
            strand = "-" if frame < 0 else "+"
    if not coords:
        print(f"{name}: NO HITS, SKIPPING")
        continue
    all_pos = [p for pair in coords for p in pair]
    lo, hi = min(all_pos), max(all_pos)
    rec = next(SeqIO.parse(os.path.join(genome_dir, f"{name}.fasta"), "fasta"))
    seq = rec.seq[lo-1:hi]
    if strand == "-":
        seq = seq.reverse_complement()
    results_nt.append((name, str(seq)))
    protein = str(Seq(str(seq)).translate(table=11, to_stop=False))
    results_aa.append((name, protein))
    print(f"{name}: {lo}-{hi} ({strand}), {len(seq)} bp -> {len(protein)} aa")
    print(f"  ** MANUALLY VERIFY unusual lengths against the genome's own GenBank annotation before trusting this extraction **")

with open("terminase_all_nt.fasta", "w") as f:
    for name, seq in results_nt: f.write(f">{name}\n{seq}\n")
with open("terminase_all_aa.fasta", "w") as f:
    for name, prot in results_aa: f.write(f">{name}\n{prot}\n")
print(f"Wrote {len(results_nt)} sequences")
PYEOF

echo "== Aligning with MAFFT =="
mafft --auto terminase_all_aa.fasta > terminase_aligned.fasta

echo "== Building tree with IQ-TREE =="
iqtree -s terminase_aligned.fasta -m MFP -bb 1000 -nt AUTO -pre terminase_tree

echo "Done. Tree: terminase_tree.treefile / terminase_tree.contree"
echo "IMPORTANT: this pipeline can mis-extract genes when a genome's terminase"
echo "diverges enough to produce fragmented tblastn hits (observed once with phiS"
echo "in the original analysis). Always cross-check unusual lengths against each"
echo "genome's own GenBank CDS annotation before trusting the extraction."
