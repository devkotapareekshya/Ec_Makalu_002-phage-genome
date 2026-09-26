#!/bin/bash
# Whole-genome comparative annotation QC: Ec_Makalu_002 vs reference MN709127.1
# Usage: run from a directory containing MN709127.gb and the Bakta annotation folder
set -e

REF_GB="MN709127.gb"
ASSEMBLY_FAA="bakta_annotation_final/Ec_Makalu_002.faa"
ASSEMBLY_FASTA="Ec_Makalu_002_final_circular.fasta"

echo "== Extracting reference proteins from GenBank =="
python3 - <<PYEOF
from Bio import SeqIO
records = list(SeqIO.parse("$REF_GB", "genbank"))
with open("MN709127_proteins.faa", "w") as out:
    count = 0
    for rec in records:
        for feature in rec.features:
            if feature.type == "CDS" and "translation" in feature.qualifiers:
                locus = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                protein_id = feature.qualifiers.get("protein_id", ["NA"])[0]
                product = feature.qualifiers.get("product", ["unknown"])[0]
                translation = feature.qualifiers["translation"][0]
                out.write(f">{locus}|{protein_id}|{product}\n{translation}\n")
                count += 1
    print(f"Extracted {count} reference proteins")
PYEOF

echo "== Building BLAST databases =="
makeblastdb -in MN709127_proteins.faa -dbtype prot -out db_reference
makeblastdb -in "$ASSEMBLY_FAA" -dbtype prot -out db_assembly

echo "== Running reciprocal BLASTp =="
blastp -query "$ASSEMBLY_FAA" -db db_reference \
    -out assembly_vs_ref_blastp.tsv \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
    -evalue 1e-5 -max_target_seqs 1 -num_threads 4

blastp -query MN709127_proteins.faa -db db_assembly \
    -out ref_vs_assembly_blastp.tsv \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
    -evalue 1e-5 -max_target_seqs 1 -num_threads 4

echo "== Identifying genes with no reciprocal hit =="
python3 - <<PYEOF
assembly_hits, ref_hits = set(), set()
with open("assembly_vs_ref_blastp.tsv") as f:
    for line in f: assembly_hits.add(line.split("\t")[0])
with open("ref_vs_assembly_blastp.tsv") as f:
    for line in f: ref_hits.add(line.split("\t")[0])

assembly_all, ref_all = set(), set()
with open("$ASSEMBLY_FAA") as f:
    for line in f:
        if line.startswith(">"): assembly_all.add(line[1:].split()[0])
with open("MN709127_proteins.faa") as f:
    for line in f:
        if line.startswith(">"): ref_all.add(line[1:].split()[0])

print("Assembly genes with no reference hit:", sorted(assembly_all - assembly_hits))
print("Reference genes with no assembly hit:", sorted(ref_all - ref_hits))
PYEOF

echo "Done. See assembly_vs_ref_blastp.tsv and ref_vs_assembly_blastp.tsv for full results."
echo "NOTE: any genes flagged above should be followed up with nucleotide-level BLASTn"
echo "against the counterpart genome before concluding a real difference exists —"
echo "see project report for the nt-level confirmation method used."
