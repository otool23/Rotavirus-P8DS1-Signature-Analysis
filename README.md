# Rotavirus P[8]DS-1 Signature Analysis

This repository contains reproducible Python scripts used to analyze genetic
signatures, allele-level genetic differentiation, and genome-wide relationships
among emergent P[8]DS-1-like Rotavirus A reassortant strains.

The repository accompanies the manuscript:

**Genetic Signatures of Emergent DS-1-like Reassortant Human Rotaviruses
Exhibiting Wa-like P[8] VP4 Genes**

The analyses include:

1. Amino-acid and nucleotide signature discovery
2. Fisher's exact testing of identified signature sites
3. Allele-level genetic-distance analysis
4. Concatenated genome alignment for evaluation of Clades I–V

---

# Repository contents

## 1. Amino-acid and nucleotide signature discovery and Fisher's exact testing

Script:

`P1_signature_discovery_with_Fisher.py`

This script identifies amino-acid and nucleotide changes associated with
reassortant-related allele groups across the analyzed Rotavirus A genome
segments.

The analysis includes:

- amino-acid signature identification
- nucleotide signature identification
- group-specific prevalence calculations
- explicit missing-data handling
- callable-site filtering
- candidate-site auditing
- Fisher's exact testing of accepted signature sites
- Benjamini-Hochberg false-discovery-rate correction
- prevalence and callable-coverage heatmaps
- Fisher-test significance heatmaps
- reproducibility metadata
- SHA-256 checksums for input FASTA files

### Genome segments analyzed

The following segments are included:

- VP1
- VP2
- VP3
- VP4
- VP6
- NSP1
- NSP2
- NSP3
- NSP4

VP7 and NSP5 are not included in these analyses.

---

## FASTA group naming

For all segments except VP4, FASTA headers use the following prefixes:

- `RP_` = P[8]DS-1 reassortant carrying allele B
- `WTP_` = wild-type strain carrying allele B
- `R_` = allele A
- `O_` = allele C
- `G_` = allele D

For VP4:

- `LGR_` = reassortant/light-green group
- `LGWT_` = wild-type strain carrying the reassortant-like/light-green VP4
- `DG_` = wild-type/deep-green VP4 group

Example headers:

`>RP_KU059768`

`>WTP_KC443224`

`>R_AB848007`

For VP4:

`>LGR_KU059769`

`>LGWT_AB848005`

`>DG_MZ546113`

---

## Input alignments

Expected amino-acid alignment names:

`VP1_AA.fasta`

`VP2_AA.fasta`

`VP3_AA.fasta`

and equivalently for the remaining analyzed segments.

Expected nucleotide alignment names:

`VP1_NU.fasta`

`VP2_NU.fasta`

`VP3_NU.fasta`

and equivalently for the remaining analyzed segments.

The FASTA files must already represent multiple-sequence alignments.

---

## Signature criteria

### Criterion 1 — strict B/light-green signature

For non-VP4 segments:

- reassortant group (`RP`) = 100% among callable sequences
- wild-type B group (`WTP`) = 100%
- allele A (`R`) = 0%
- allele C (`O`) = 0%, if present
- allele D (`G`) = 0%, if present

For VP4:

- `LGR` = 100%
- `LGWT` = 100%
- `DG` = 0%

### Criterion 2 — reassortant-enriched B/light-green signature

For non-VP4 segments:

- `RP` = 100%
- `WTP` > 0% but < 100%
- `R`, `O`, and `G` = 0%, when those groups are present

For VP4:

- `LGR` = 100%
- `LGWT` > 0% but < 100%
- `DG` = 0%

By default, an existing required group must have at least 80% callable
sequences at a site before the site can satisfy a signature criterion.

Groups that are completely absent for a gene are treated as not applicable
rather than as evidence against a signature.

---

## Missing-data handling

Gaps and ambiguous amino-acid or nucleotide states are treated as
non-callable.

Signature prevalence is calculated using:

`State prevalence (%) = n_state / n_callable × 100`

where:

- `n_state` = number of callable sequences carrying the candidate state
- `n_callable` = number of sequences with an unambiguous amino acid or
  nucleotide at that position

Missing or ambiguous calls are therefore excluded from the prevalence
denominator but are reported separately.

---

## Mutation notation

Mutations are reported as:

`WT state + biological position + reassortant state`

Examples:

`V296I`

`Y457C`

`A742G`

The biological position is based on the ungapped wild-type reference sequence
selected by the script.

---

# Fisher's exact tests

Fisher's exact tests are performed for every accepted Criterion 1 and
Criterion 2 signature site.

Only callable sequences are included in the contingency tables.

The general 2 × 2 table is:

| Group | Signature state | Other callable state |
|---|---:|---:|
| Reassortant/test group | a | b |
| Comparison group | c | d |

A two-sided Fisher's exact test is performed.

## Non-VP4 comparisons

Two types of tests are performed:

### Reassortant-only analysis

- `RP` vs `R`
- `RP` vs `O`, when present
- `RP` vs `G`, when present

### Combined B-associated analysis

- `RP + WTP` vs `R`
- `RP + WTP` vs `O`, when present
- `RP + WTP` vs `G`, when present

## VP4 comparisons

- `LGR` vs `DG`
- `LGR + LGWT` vs `DG`

These tests are performed for both Criterion 1 and Criterion 2 sites.

The outputs include:

- explicit numerators and denominators
- odds ratio
- two-sided Fisher's exact p-value
- Benjamini-Hochberg FDR-adjusted p-value

Benjamini-Hochberg correction is reported both within predefined comparison
families and globally across all Fisher tests.

The within-family FDR is the primary corrected value used for interpretation.

### Fisher-test outputs

`P1_signature_Fisher_tests.csv`

`P1_AA_signature_Fisher_tests.csv`

`P1_NU_signature_Fisher_tests.csv`

Fisher-test heatmaps are also generated separately for amino-acid and
nucleotide signatures.

---

# 2. Allele genetic-distance analysis

Script:

`P1_allele_genetic_distance_analysis.py`

This analysis quantitatively evaluates genetic differentiation between allele B
and the alternative intra-genotypic allele groups.

Expected pairwise percent-identity matrices are:

`VP1.csv`

`VP2.csv`

`VP3.csv`

`VP4.csv`

`VP6.csv`

`NSP1.csv`

`NSP2.csv`

`NSP3.csv`

`NSP4.csv`

Sequence labels in these matrices begin with:

- `A-` = allele A
- `B-` = allele B
- `C-` = allele C
- `D-` = allele D

Example:

`A-KJ123456`

`B-KC443369`

Allele B is compared with alleles A, C, and D when those groups are present.

---

## Genetic-distance calculations

Pairwise genetic distance is calculated as:

`p-distance (%) = 100 - percent identity`

For each comparison, the script calculates:

- mean within-B identity
- mean within-other-allele identity
- mean between-group identity
- mean within-B genetic distance
- mean within-other-allele genetic distance
- mean between-group genetic distance
- net between-group genetic distance
- identity ratio
- distance ratio

Net between-group distance is calculated as:

`d_net = d_between - (d_within_B + d_within_other) / 2`

Identity ratio:

`within-B identity / between-group identity`

Distance ratio:

`within-B distance / between-group distance`

An identity ratio greater than 1 indicates that allele B sequences are more
similar to one another than they are to the comparison allele.

A distance ratio below 1 indicates that within-B genetic divergence is smaller
than the divergence between allele B and the comparison allele.

These quantitative analyses were used to support the intra-genotypic allele
assignments initially observed in the segment-specific phylogenetic trees.

---

# 3. Concatenated genome analysis for evaluation of Clades I–V

Script:

`P1_build_concatenated_reassortant_MSA_NCBI_fallback.py`

A concatenated nucleotide alignment was generated to evaluate the genome-wide
relationships among reassortant strains assigned to Clades I–V.

For each reassortant strain, the following segments are combined in a fixed
order:

`VP4 | VP6 | VP1 | VP2 | VP3 | NSP1 | NSP2 | NSP3 | NSP4`

Each segment is aligned separately using MAFFT before concatenation.

The same set of strains is retained for every segment, ensuring that each row
of the final concatenated alignment represents one reassortant strain.

The final output is:

`P1_reassortant_concatenated_MSA.fasta`

Each FASTA record therefore represents one concatenated reassortant genome
across the nine analyzed segments.

---

## Accession table

The concatenation analysis uses a CSV table containing one row per strain and
one accession number per segment.

Required columns are:

`Clade,Strain,VP4,VP6,VP1,VP2,VP3,NSP1,NSP2,NSP3,NSP4`

The accession table used for this analysis is:

`Clades_accessions.csv`

---

## Automatic NCBI fallback

If an accession listed in the clade table is not found in the corresponding
local nucleotide FASTA file, the script retrieves the exact nucleotide record
from NCBI.

Downloaded sequences are saved separately and the original local FASTA files
are not modified.

The downloaded sequence is then included with the corresponding segment
sequences before MAFFT alignment and subsequent concatenation.

---

## Clade evaluation

The concatenated whole-genome analysis was used as an independent evaluation
of the previously defined constellation-based Clades I–V.

The resulting genome-wide phylogenetic analysis generally supported the
original clade assignments, with strains within each clade clustering together.

Three strains did not cluster with their initially assigned groups and were
treated as exceptions because of their atypical allele constellations:

- `IRL/16IRL33182/2016/G3P[8]`
- `THA/DBM2018-291/2018/G9P[8]`
- `DOM/3000503734/2016/G3P[8]`

---

# Requirements

Python 3

Required Python packages:

- `numpy`
- `matplotlib`
- `scipy`
- `biopython`

Install the Python dependencies using:

`pip install -r requirements.txt`

MAFFT is additionally required for the concatenated genome analysis and must
be installed separately or available through the computing environment.

The concatenation workflow was run using MAFFT 7.525.

---

# Running the analyses

## Signature discovery and Fisher tests

`python P1_signature_discovery_with_Fisher.py`

Optional:

`python P1_signature_discovery_with_Fisher.py --input-dir . --output-dir P1_signature_results --min-callable-fraction 0.80`

---

## Allele genetic-distance analysis

`python P1_allele_genetic_distance_analysis.py`

Optional:

`python P1_allele_genetic_distance_analysis.py --input-dir . --output-dir P1_distance_results`

---

## Concatenated reassortant genome alignment

Example:

`python P1_build_concatenated_reassortant_MSA_NCBI_fallback.py --accession-table Clades_accessions.csv --input-dir . --output-dir concatenated_MSA --ncbi-email your_email@example.com`

MAFFT must be available in the system PATH before running this analysis.

---

# Reproducibility

The analyses were designed to facilitate independent reproduction of the
results.

Where applicable, the scripts record:

- software versions
- analysis parameters
- sequence group counts
- input-file checksums
- callable and missing sequence counts
- accepted signature sites
- rejected candidate sites
- Fisher-test contingency-table counts
- raw and FDR-adjusted p-values
- sequence sources used for concatenation
- segment boundaries in the concatenated alignment

No random sampling is used in the signature-discovery or genetic-distance
analyses.

---

# Data and code availability

The analysis scripts used in this study are publicly available at:

https://github.com/otool23/Rotavirus-P8DS1-Signature-Analysis

---

# Citation

If you use this code, please cite the associated manuscript:

**Genetic Signatures of Emergent DS-1-like Reassortant Human Rotaviruses
Exhibiting Wa-like P[8] VP4 Genes**

Full citation will be added following publication.

---

# License

This repository is distributed under the MIT License.
