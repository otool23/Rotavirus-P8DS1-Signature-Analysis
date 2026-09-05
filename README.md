# Rotavirus P[8]DS-1 Signature Analysis

This repository contains reproducible Python scripts used in the analysis of
genetic signatures and allele-level genetic differentiation in emergent
P[8]DS-1-like Rotavirus A reassortant strains.

The repository accompanies the manuscript:

**Genetic Signatures of Emergent DS-1-like Reassortant Human Rotaviruses
Exhibiting Wa-like P[8] VP4 Genes**

## Repository contents

### 1. Amino-acid and nucleotide signature discovery

Script:

`P1_signature_discovery_reproducible.py`

This script identifies amino-acid and nucleotide changes associated with
reassortant-related allele groups across the analyzed Rotavirus A genome
segments.

The analysis includes:

- amino-acid signature identification
- nucleotide signature identification
- group-specific prevalence calculations
- missing-data handling
- callable-site filtering
- candidate-site auditing
- generation of prevalence and callable-coverage heatmaps
- reproducibility metadata and input-file checksums

For most genome segments, sequence headers are assigned using the following
prefixes:

- `RP_` = reassortant strains carrying the B/pink allele
- `WTP_` = wild-type strains carrying the B/pink allele
- `R_` = type A/red allele
- `O_` = type C/orange allele
- `G_` = type D/grey allele

For VP4:

- `LGR_` = reassortant strains carrying the light-green VP4 allele
- `LGWT_` = wild-type strains carrying the light-green VP4 allele
- `DG_` = wild-type/deep-green VP4 allele

Expected amino-acid alignment names include:

`VP1_AA.fasta`, `VP2_AA.fasta`, etc.

Expected nucleotide alignment names include:

`VP1_NU.fasta`, `VP2_NU.fasta`, etc.

VP7 and NSP5 are not included in this analysis.

### Signature criteria

Criterion 1 identifies sites where the reassortant-associated state is fixed
among callable sequences in both the reassortant and corresponding wild-type
B-like groups and absent from the comparison allele groups.

Criterion 2 identifies sites where the reassortant-associated state is fixed
in the reassortant group but variable in the corresponding wild-type B-like
group and absent from the comparison allele groups.

By default, an existing group must have at least 80% callable sequences at a
site before the site can satisfy the signature criteria.

Missing, gapped, or ambiguous states are excluded from the prevalence
denominator.

---

### 2. Allele genetic-distance analysis

Script:

`P1_allele_genetic_distance_analysis.py`

This script quantifies genetic differentiation between allele B and the other
allele groups using pairwise percent-identity matrices.

Expected input files are:

`VP1.csv`
`VP2.csv`
`VP3.csv`
`VP4.csv`
`VP6.csv`
`NSP1.csv`
`NSP2.csv`
`NSP3.csv`
`NSP4.csv`

Sequence labels in each identity matrix must begin with:

- `A-` = allele A
- `B-` = allele B
- `C-` = allele C
- `D-` = allele D

For example:

`A-KJ123456`

`B-KC443369`

The script compares allele B with allele A, C, and D when those groups are
present.

### Genetic-distance calculations

Pairwise genetic distance is calculated as:

`p-distance (%) = 100 - percent identity`

The following statistics are calculated for each allele comparison:

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

Identity ratio is calculated as:

`within-B identity / between-group identity`

Distance ratio is calculated as:

`within-B distance / between-group distance`

An identity ratio greater than 1 indicates that allele B sequences are more
similar to one another than they are to the comparison allele.

A distance ratio below 1 indicates that within-B genetic divergence is smaller
than the divergence between allele B and the comparison allele.

---

## Requirements

Python 3

Required Python packages:

`numpy`

`matplotlib`

Install the dependencies using:

`pip install -r requirements.txt`

---

## Running the analyses

### Signature analysis

`python P1_signature_discovery_reproducible.py`

### Allele genetic-distance analysis

`python P1_allele_genetic_distance_analysis.py`

The scripts may also be run with explicit input and output directories.

Example:

`python P1_allele_genetic_distance_analysis.py --input-dir . --output-dir P1_distance_results`

---

## Reproducibility

Both scripts were designed to support reproducible analysis.

The scripts record information such as:

- software versions
- analysis parameters
- input-file SHA-256 checksums
- output files
- group counts

No random sampling is used in the allele genetic-distance analysis.

---

## Citation

If you use this code, please cite the associated manuscript:

[Insert full manuscript citation after publication]

---

## License

This repository is distributed under the MIT License.
