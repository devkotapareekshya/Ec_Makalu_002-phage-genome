"""
Extract CDS protein translations from a GenBank record into a FASTA file,
for use as the reference proteome in reciprocal BLASTp comparative QC.

Injected by Snakemake: snakemake.input.gb, snakemake.output.faa
"""

from Bio import SeqIO

gb_path = snakemake.input.gb
faa_path = snakemake.output.faa

n_written = 0
with open(faa_path, "w") as out:
    for record in SeqIO.parse(gb_path, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            translation = feature.qualifiers.get("translation")
            if not translation:
                continue
            locus_tag = feature.qualifiers.get("locus_tag", ["unknown"])[0]
            product = feature.qualifiers.get("product", [""])[0]
            out.write(f">{locus_tag} {product}\n{translation[0]}\n")
            n_written += 1

print(f"Extracted {n_written} protein translations from {gb_path} -> {faa_path}")
