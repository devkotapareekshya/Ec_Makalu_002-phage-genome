"""
Circularize the longest SPAdes contig and reorient it so its start coordinate
matches the reference genome's start (for easier comparative annotation).

This is a Snakemake `script:` directive -- `snakemake` object is injected
automatically with .input, .output, .params, .threads.

Approach:
  1. Take the longest contig from the SPAdes assembly.
  2. Check for terminal overlap between the contig's start and end (a common
     signature of a circular genome that SPAdes assembled as a linear contig
     with redundant terminal sequence). If found, trim the redundant copy.
  3. BLAST the (now-circularized) contig against the reference to find the
     reference's start position on the query, then rotate the sequence so
     that position becomes position 1.

CAUTION: this heuristic works well for phage genomes with a single dominant
contig and a close reference, which is the case here, but always sanity-check
the output (e.g. re-run the terminal_repeat_check rule) before trusting it
blindly on a new genome.
"""

import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

contigs_fasta = snakemake.input.contigs
ref_fasta = snakemake.input.ref
out_fasta = snakemake.output.final

workdir = Path(out_fasta).parent
workdir.mkdir(parents=True, exist_ok=True)

# --- 1. Pick the longest contig ---------------------------------------------
records = list(SeqIO.parse(contigs_fasta, "fasta"))
if not records:
    sys.exit(f"No contigs found in {contigs_fasta}")
longest = max(records, key=lambda r: len(r.seq))
seq = str(longest.seq).upper()
print(f"Longest contig: {longest.id}, length {len(seq)} bp")

# --- 2. Check for terminal redundancy (naive exact-match overlap scan) -----
# Look for the largest exact match between a suffix and a prefix of the contig.
max_check = min(2000, len(seq) // 4)
overlap = 0
for k in range(max_check, 10, -1):
    if seq[:k] == seq[-k:]:
        overlap = k
        break

if overlap > 0:
    print(f"Detected {overlap} bp terminal redundancy -- trimming duplicate copy.")
    seq = seq[:-overlap]
else:
    print("No terminal redundancy detected via exact-match scan; "
          "assuming contig is already a clean single-copy representation.")

circ_fasta = workdir / "circularized_untrimmed.fasta"
SeqIO.write(
    SeqRecord(Seq(seq), id=f"{longest.id}_circular", description=""),
    circ_fasta, "fasta",
)

# --- 3. Reorient to match reference start -----------------------------------
# BLAST the circularized contig against the reference; find where the
# reference's position 1 maps onto our contig, then rotate.
db_prefix = str(workdir / "circ_db")
subprocess.run(
    ["makeblastdb", "-in", str(circ_fasta), "-dbtype", "nucl", "-out", db_prefix],
    check=True,
)

# Use the first 200 bp of the reference as an anchor query.
ref_records = list(SeqIO.parse(ref_fasta, "fasta"))
ref_start_anchor = workdir / "ref_start_anchor.fasta"
SeqIO.write(
    SeqRecord(ref_records[0].seq[:200], id="ref_anchor", description=""),
    ref_start_anchor, "fasta",
)

blast_out = workdir / "anchor_hit.tsv"
subprocess.run(
    [
        "blastn", "-query", str(ref_start_anchor), "-db", db_prefix,
        "-outfmt", "6 sstart send sstrand pident length",
        "-max_target_seqs", "1", "-out", str(blast_out),
    ],
    check=True,
)

rotate_pos = 0
with open(blast_out) as fh:
    line = fh.readline().strip()
    if line:
        sstart, send, sstrand, pident, length = line.split("\t")
        rotate_pos = int(sstart) - 1  # 0-based rotation point
        print(f"Reference start anchor maps to contig position {sstart} "
              f"(strand={sstrand}, identity={pident}%, length={length})")
    else:
        print("WARNING: no BLAST hit for reference start anchor -- "
              "leaving contig unrotated. Inspect manually.")

if rotate_pos > 0:
    seq = seq[rotate_pos:] + seq[:rotate_pos]

final_id = snakemake.wildcards.get("sample", "assembly") if hasattr(snakemake, "wildcards") else "assembly"
final_record = SeqRecord(Seq(seq), id=f"{final_id}_final_circular", description="")
SeqIO.write(final_record, out_fasta, "fasta")
print(f"Wrote final circularized, reoriented genome to {out_fasta} "
      f"({len(seq)} bp)")
