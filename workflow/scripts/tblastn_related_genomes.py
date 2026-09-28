"""
Run tblastn of the Ec_Makalu_002 TerL protein against each related-taxon
genome FASTA in the phylogenetics/phylo_genomes/ directory, extract the
best-hit region, translate it, and collect all homologs into one FASTA
ready for MAFFT alignment.

Injected by Snakemake: snakemake.input.query, snakemake.input.genomes_dir,
snakemake.output.hits, snakemake.threads
"""

import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

query = snakemake.input.query
genomes_dir = Path(snakemake.input.genomes_dir)
out_faa = snakemake.output.hits
threads = snakemake.threads

genome_files = sorted(genomes_dir.glob("*.fasta")) + sorted(genomes_dir.glob("*.fa"))
if not genome_files:
    raise FileNotFoundError(
        f"No FASTA files found in {genomes_dir}. Populate this directory "
        "with related Krischvirus/Straboviridae genome FASTA files before "
        "running the phylogenetics rules (see report Methods 2.3 for the "
        "15 taxa used in the original analysis)."
    )

records_out = []

for genome_fasta in genome_files:
    taxon = genome_fasta.stem
    db_prefix = f"/tmp/{taxon}_db"
    subprocess.run(
        ["makeblastdb", "-in", str(genome_fasta), "-dbtype", "nucl", "-out", db_prefix],
        check=True, capture_output=True,
    )

    result = subprocess.run(
        [
            "tblastn", "-query", query, "-db", db_prefix,
            "-outfmt", "6 sseqid sstart send sstrand evalue bitscore",
            "-max_target_seqs", "1", "-num_threads", str(threads),
        ],
        check=True, capture_output=True, text=True,
    )

    line = result.stdout.strip().split("\n")[0] if result.stdout.strip() else None
    if not line:
        print(f"WARNING: no tblastn hit for {taxon} -- skipping.")
        continue

    sseqid, sstart, send, sstrand, evalue, bitscore = line.split("\t")
    sstart, send = int(sstart), int(send)
    low, high = min(sstart, send), max(sstart, send)

    genome_records = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    seq = genome_records[sseqid].seq[low - 1:high]
    if sstrand == "minus":
        seq = seq.reverse_complement()

    protein = seq.translate(table=11, to_stop=True)
    records_out.append(SeqRecord(protein, id=taxon, description=f"tblastn_hit_{sseqid}:{low}-{high}"))
    print(f"{taxon}: hit at {sseqid}:{low}-{high} ({sstrand}), "
          f"evalue={evalue}, bitscore={bitscore}, translated length={len(protein)} aa")

SeqIO.write(records_out, out_faa, "fasta")
print(f"Wrote {len(records_out)} homologous TerL sequences -> {out_faa}")
