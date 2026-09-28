#!/usr/bin/env python3
"""
Parse forward and reverse BLASTp outfmt-6 results and report reciprocal best
hits (RBH) between the query and reference proteomes, plus any proteins on
either side that lack a reciprocal hit (candidates for the nucleotide-level
BLASTn follow-up described in the report's comparative annotation QC).

Usage:
    python3 reciprocal_best_hits.py fwd.tsv rev.tsv output.tsv

fwd.tsv: query_vs_ref blastp outfmt 6 (query proteome vs. reference db)
rev.tsv: ref_vs_query blastp outfmt 6 (reference proteome vs. query db)
"""

import sys
from collections import defaultdict

fwd_path, rev_path, out_path = sys.argv[1:4]


def best_hits(path):
    """Return {query_id: best_subject_id} keeping only the top hit per query
    (outfmt 6 rows are already sorted best-first per query by blastp with
    -max_target_seqs 1, but we defensively keep the first row seen)."""
    hits = {}
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            qid, sid = fields[0], fields[1]
            if qid not in hits:
                hits[qid] = sid
    return hits


fwd = best_hits(fwd_path)   # query -> ref
rev = best_hits(rev_path)   # ref -> query

rbh = []
query_no_hit = []
ref_no_hit = []

for q, r in fwd.items():
    if rev.get(r) == q:
        rbh.append((q, r))

all_queries = set(fwd.keys())
all_refs = set(rev.keys())
rbh_queries = {q for q, r in rbh}
rbh_refs = {r for q, r in rbh}

query_no_hit = sorted(all_queries - rbh_queries)
ref_no_hit = sorted(all_refs - rbh_refs)

with open(out_path, "w") as out:
    out.write("type\tquery_id\tref_id\n")
    for q, r in rbh:
        out.write(f"RBH\t{q}\t{r}\n")
    for q in query_no_hit:
        out.write(f"QUERY_NO_RECIPROCAL_HIT\t{q}\t{fwd.get(q,'')}\n")
    for r in ref_no_hit:
        out.write(f"REF_NO_RECIPROCAL_HIT\t{rev.get(r,'')}\t{r}\n")

print(f"Reciprocal best hits: {len(rbh)}")
print(f"Query proteins with no reciprocal hit: {len(query_no_hit)}")
print(f"Reference proteins with no reciprocal hit: {len(ref_no_hit)}")
print("NOTE: proteins with no reciprocal hit should be followed up with "
      "nucleotide-level BLASTn against the counterpart genome before "
      "concluding they represent genuine gene-content differences -- see "
      "report Section 3.4 for why this matters (annotation-pipeline "
      "sensitivity artifacts vs. real divergence).")
