# *Escherichia* phage Ec_Makalu_002 — Genome Annotation and Comparative Analysis

Complete genome assembly, annotation, and comparative/phylogenetic analysis of a novel *Escherichia* phage in the genus *Krischvirus* (family *Straboviridae*, order *Pantevenvirales* — the modern classification containing T4 and its relatives).

## Summary

- **164,673 bp** circular genome, essentially identical in length to the closest public reference (MN709127.1, 164,674 bp)
- **269 predicted CDS**, annotated with Bakta
- Whole-genome comparison against the reference confirmed **no genuine gene-content differences** — all apparent discrepancies traced to differences in ORF-calling sensitivity between annotation pipelines, not real biological divergence
- Phylogenetic analysis of the terminase large subunit placed Ec_Makalu_002 within a distinct subclade of *Krischvirus*
- That analysis led to a novel finding: a **~330 amino acid insertion** in the terminase large subunit, shared by this subclade and bearing sequence hallmarks of a **mobile group I intein** (conserved N-terminal splice residue, canonical C-terminal His-Asn splicing motif, and a histidine cluster consistent with an embedded HNH homing endonuclease domain)

See [`report/Ec_Makalu_002_mini_report.docx`](report/Ec_Makalu_002_mini_report.docx) for the full write-up with methods, figures, and discussion.

## Key findings

### 1. Long tail fiber gene (host-range determinant)
The long tail fiber gene was fully resolved at its correct coordinates and shows 99.92% identity to the reference (1246/1246 aa aligned, single conservative substitution at position 245).

### 2. Whole-genome comparative annotation QC
Reciprocal BLASTp against the reference genome's annotated proteome (274 CDS) found 268/269 assembly genes and 268/274 reference genes with confident reciprocal hits. The 7 remaining apparent discrepancies were each resolved by nucleotide-level BLAST: **all 7 showed 100% identity, full-length, zero mismatches** against the counterpart genome. Conclusion: the underlying DNA is fully conserved; the discrepancies are gene-calling threshold artifacts, not true gaps or novel genes.

### 3. Putative mobile intein in the terminase large subunit
A phylogenetic analysis of the terminase large subunit gene across 16 *Krischvirus* and related genomes (MAFFT alignment + IQ-TREE, ultrafast bootstrap) showed Ec_Makalu_002 clustering with Ec_Makalu_003, ECD7, and AV108 in a well-supported subclade (bootstrap ≥ 89) distinguished by a ~330 aa insertion absent from the rest of the genus. This insertion:
- Begins with Ala, matching the canonical intein N-terminal splice-junction residue
- Ends in the conserved His-Asn dipeptide (Block G), the canonical intein C-terminal splicing motif
- Contains a histidine cluster in its central third, consistent with an embedded HNH-family homing endonuclease domain
- Shares its insertion site (221/222 overlapping alignment columns) with a shorter, independent 222 aa insertion found in the more distantly related phage suisiam, suggesting this is a recurrent intein insertion hotspot in the *Krischvirus* terminase large subunit gene

This finding is presented as a well-supported sequence-based hypothesis; formal confirmation (e.g., against InBase or a dedicated intein-splicing domain profile) was not completed due to intermittent NCBI service availability, and is noted as a next step.

### 4. GLGIAI_0263 — conservatively labeled
 true gaps or novel genes.

### 3. Putative mobile intein in the terminase large subunit
A phylogenetic analysis of the terminase large subunit gene across 16 *Krischvirus* and related genomes (MAFFT alignment + IQ-TREE, ultrafast bootstrap) showed Ec_Makalu_002 clustering with Ec_Makalu_003, ECD7, and AV108 in a well-supported subclade (bootstrap ≥ 89) distinguished by a ~330 aa insertion absent from the rest of the genus. This insertion:
- Begins with Ala, matching the canonical intein N-terminal splice-junction residue
- Ends in the conserved His-Asn dipeptide (Block G), the canonical intein C-terminal splicing motif
- Contains a histidine cluster in its central third, consistent with an embedded HNH-family homing endonuclease domain
- Shares its insertion site (221/222 overlapping alignment columns) with a shorter, independent 222 aa insertion found in the more distantly related phage suisiam, suggesting this is a recurrent intein insertion hotspot in the *Krischvirus* terminase large subunit gene

This finding is presented as a well-supported sequence-based hypothesis; formal confirmation (e.g., against InBase or a dedicated intein-splicing domain profile) was not completed due to intermittent NCBI service availability, and is noted as a next step.

### 4. GLGIAI_0263 — conservatively labeled
A 109 aa CDS with only distant homology to an uncharacterized ECD7 protein and no confident Pfam-A domain hit is labeled `putative homing endonuclease` throughout the annotation, consistent with the original NCBI submission's own qualified call for the same locus.

## Repository structure

├── data/ Final circularized assembly (164,673 bp)
├── annotation/bakta_annotation_final/ Full Bakta annotation output
├── comparative_analysis/
│ ├── scripts/run_comparison.sh Reciprocal BLASTp pipeline
│ └── results/ BLAST tables
├── phylogenetics/
│ ├── scripts/build_tree.sh Extraction + alignment + tree pipeline
│ ├── phylo_genomes/ 16 taxon genomes used
│ └── results/ Alignment, tree files
└── report/
└── Ec_Makalu_002_mini_report.docx Full write-up


## Methods overview

- **Annotation**: Bakta, default parameters
- **Comparative genomics**: reciprocal BLASTp (E-value ≤ 1e-5) between Bakta CDS calls and reference GenBank CDS translations; discrepancies resolved with nucleotide-level BLASTn
- **Phylogenetics**: terminase large subunit identified via Pfam domains (PF14528, PF17289); homologs recovered from 15 related genomes via tblastn; aligned with MAFFT; tree built with IQ-TREE (ModelFinder + 1000 ultrafast bootstrap replicates)
- **Reference genome**: GenBank accession MN709127.1

## Reproducing this analysis

```bash
cd comparative_analysis/scripts && bash run_comparison.sh
cd phylogenetics/scripts && bash build_tree.sh
```

## Notes on provenance

This project was developed iteratively, including using AI-assisted analysis for pipeline execution, comparative BLAST interpretation, and phylogenetic tree construction. All scientific conclusions were verified against primary evidence rather than taken at face value, and errors identified during the process (e.g., an incorrect tail fiber gene coordinate from an earlier working session) were corrected before finalizing results.
