# Rotavirus A P[8]DS-1 Reassortant Signature Analysis

This repository contains the reproducible Python analysis used to identify
amino-acid and nucleotide signatures associated with emergent P[8]DS-1-like
Rotavirus A reassortant strains.

## Requirements

Python 3.10 or later

Python packages:
- numpy
- matplotlib

## Input files

The script analyzes aligned amino-acid and nucleotide FASTA files:

VP1_AA.fasta
VP1_NU.fasta
VP2_AA.fasta
VP2_NU.fasta
...

VP7 and NSP5 are not included in this analysis.

## FASTA header naming

For VP1-VP6 and NSP1-NSP4:

RP_   = P[8]DS-1 reassortant type-B allele
WTP_  = wild-type strain carrying type-B allele
R_    = type-A allele
O_    = type-C allele
G_    = type-D allele

Examples:

>RP_KJ123456
>WTP_KC443224
>R_KC443369

For VP4:

LGR_  = reassortant/light-green group
LGWT_ = wild-type light-green group
DG_   = wild-type/deep-green group

## Signature criteria

Criterion 1:
100% of callable RP sequences and 100% of callable WTP sequences carry
the reassortant state, while the state is absent from callable R, O, and G
sequences.

Criterion 2:
100% of callable RP sequences carry the reassortant state, the state is
present but variable in WTP sequences, and absent from callable R, O, and G
sequences.

For VP4, the corresponding groups are LGR, LGWT, and DG.

## Missing data

Gaps and ambiguous residues/bases are treated as missing data and excluded
from site-specific denominators. The script reports total N, callable N,
missing N, and callable fraction for each group.

The default minimum callable fraction is 0.80.

## Running the analysis

python P1_signature_discovery_reproducible.py

or:

python P1_signature_discovery_reproducible.py \
    --input-dir . \
    --output-dir P1_signature_results \
    --min-callable-fraction 0.80

## Outputs

The analysis generates:

- P1_AA_signature_changes.csv
- P1_NU_signature_changes.csv
- P1_signature_candidate_audit.csv
- P1_signature_failure_details.csv
- P1_signature_group_counts.csv
- P1_signature_run_metadata.txt
- prevalence and callable-coverage figures
