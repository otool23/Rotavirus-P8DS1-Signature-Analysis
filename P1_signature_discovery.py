#!/usr/bin/env python3
"""
Paper 1: reproducible reassortant-signature discovery
=====================================================

Purpose
-------
Identify amino-acid (AA) and nucleotide (NU) signatures associated with the
P[8]DS-1 reassortant / type-B background using predefined FASTA-header groups.

This version explicitly handles missing sequence calls (gaps and ambiguous
characters) at each aligned site. Missing calls are EXCLUDED from the
site-specific percentage denominator, but are reported separately as total N,
callable N, missing N, and callable fraction.

Expected input files
--------------------
Placed in the input directory (default: directory containing this script):

    VP1_AA.fasta   VP1_NU.fasta
    VP2_AA.fasta   VP2_NU.fasta
    VP3_AA.fasta   VP3_NU.fasta
    VP4_AA.fasta   VP4_NU.fasta
    VP6_AA.fasta   VP6_NU.fasta
    NSP1_AA.fasta  NSP1_NU.fasta
    NSP2_AA.fasta  NSP2_NU.fasta
    NSP3_AA.fasta  NSP3_NU.fasta
    NSP4_AA.fasta  NSP4_NU.fasta

VP7 and NSP5 are intentionally excluded.

Header-prefix groups
--------------------
Non-VP4:
    RP_   P[8]DS-1 reassortant, type-B/pink allele
    WTP_  wild-type strain carrying type-B/pink allele
    R_    wild-type type-A/red allele
    O_    type-C/orange allele
    G_    type-D/grey allele

VP4:
    LGR_  reassortant/light-green group
    LGWT_ wild-type strain carrying light-green/reassortant-like VP4
    DG_   wild-type/deep-green group

Signature criteria
------------------
Criterion 1 -- strict B/light-green signature

    Non-VP4:
        RP  = 100% among callable sequences
        WTP = 100% among callable sequences
        R/O/G = 0% among callable sequences, IF that group exists

    VP4:
        LGR  = 100% among callable sequences
        LGWT = 100% among callable sequences
        DG   = 0% among callable sequences

Criterion 2 -- reassortant-enriched B/light-green signature

    Non-VP4:
        RP  = 100% among callable sequences
        0% < WTP < 100% among callable sequences
        R/O/G = 0% among callable sequences, IF that group exists

    VP4:
        LGR  = 100% among callable sequences
        0% < LGWT < 100% among callable sequences
        DG   = 0% among callable sequences

Missing-data rule
-----------------
At a site, '-', '?', and ambiguous symbols are treated as missing/non-callable.
Percentages are calculated as:

    percent = state_count / callable_N * 100

NOT state_count / total_group_N.

To avoid calling a signature from very sparse data, a group that is required
for a criterion must meet --min-callable-fraction (default 0.80 = 80%).
A group that does not exist at all (e.g., no G_ sequences for a gene) is NA and
does not invalidate the signature. A group that exists but has insufficient
callable data at that site is NOT treated as evidence of 0%.

Mutation naming
---------------
Mutation labels use:

    WT_state + biological_position + reassortant_state

Examples:
    V296I
    Y457C
    A742G

WT state is the modal callable state in R_ (or DG_ for VP4).
Position numbering is based on a deterministic ungapped WT reference sequence:
the longest R_ sequence (or DG_ for VP4), breaking ties alphabetically by FASTA
header. The selected reference is reported in the output.


Fisher exact-test extension
---------------------------
For every accepted Criterion 1 and Criterion 2 signature site, this version
also performs two-sided Fisher's exact tests using callable sequences only.

Non-VP4:
    - RP only versus each existing comparison group (R, O, G)
    - RP + WTP combined versus each existing comparison group (R, O, G)

VP4:
    - LGR only versus DG
    - LGR + LGWT combined versus DG

Benjamini-Hochberg FDR correction is reported both within each predefined test
family and globally across all Fisher tests. Separate combined, AA-only, and
NU-only Fisher CSV files and AA/NU FDR heatmaps are produced.

Reproducibility outputs
-----------------------
The script records:
    - analysis parameters
    - Python / NumPy / Matplotlib versions
    - SHA-256 checksums for every input FASTA
    - exact group counts
    - accepted signatures
    - rejected candidate sites and reasons
    - sequence-level failure details
    - prevalence and coverage heatmaps

No random sampling is used; repeated runs on identical inputs and parameters
produce identical tabular results.

Example
-------
    python P1_signature_discovery_reproducible.py

Optional:
    python P1_signature_discovery_reproducible.py \
        --input-dir . \
        --output-dir results \
        --min-callable-fraction 0.80
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import platform
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import fisher_exact


# ============================================================================
# CONSTANTS
# ============================================================================

SCRIPT_VERSION = "2.1.0"

SEGMENTS = [
    "VP1", "VP2", "VP3", "VP4", "VP6",
    "NSP1", "NSP2", "NSP3", "NSP4",
]
SEQ_TYPES = ["AA", "NU"]

NON_VP4_GROUPS = ["RP", "WTP", "R", "O", "G"]
VP4_GROUPS = ["LGR", "LGWT", "DG"]

# Unambiguous biological states accepted as callable.
AA_VALID = set("ACDEFGHIKLMNPQRSTVWY")
NU_VALID = set("ACGTU")


# ============================================================================
# COMMAND-LINE ARGUMENTS
# ============================================================================

def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(
        description=(
            "Discover reassortant-associated AA and nucleotide signatures "
            "with explicit missing-data handling and reproducibility metadata."
        )
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        default=script_dir,
        help="Directory containing *_AA.fasta and *_NU.fasta files "
             "(default: script directory).",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=script_dir / "P1_signature_results",
        help="Directory for analysis outputs "
             "(default: ./P1_signature_results beside script).",
    )

    parser.add_argument(
        "--min-callable-fraction",
        type=float,
        default=0.80,
        help=(
            "Minimum fraction of sequences in an existing required group that "
            "must have an unambiguous state at a site (default: 0.80)."
        ),
    )

    args = parser.parse_args()

    if not 0.0 < args.min_callable_fraction <= 1.0:
        parser.error("--min-callable-fraction must be >0 and <=1.")

    args.input_dir = args.input_dir.resolve()
    args.output_dir = args.output_dir.resolve()

    return args


# ============================================================================
# FILE / REPRODUCIBILITY HELPERS
# ============================================================================

def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_run_metadata(
    output_dir: Path,
    input_files: list[Path],
    min_callable_fraction: float,
) -> None:
    """Write software versions, parameters, and input checksums."""
    metadata_path = output_dir / "P1_signature_run_metadata.txt"

    lines = [
        "Paper 1 signature discovery: reproducibility metadata",
        "====================================================",
        f"Script version: {SCRIPT_VERSION}",
        f"Run timestamp (UTC): {datetime.now(timezone.utc).isoformat()}",
        f"Python version: {platform.python_version()}",
        f"Python executable: {sys.executable}",
        f"Platform: {platform.platform()}",
        f"NumPy version: {np.__version__}",
        f"Matplotlib version: {matplotlib.__version__}",
        f"SciPy version: {scipy.__version__}",
        f"Minimum callable fraction: {min_callable_fraction:.4f}",
        "",
        "Input FASTA SHA-256 checksums",
        "----------------------------",
    ]

    for path in sorted(input_files, key=lambda p: p.name):
        lines.append(f"{path.name}\t{sha256_file(path)}")

    metadata_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ============================================================================
# FASTA AND GROUPING
# ============================================================================

def classify_header(header: str, segment: str) -> str | None:
    """Assign a sequence to its prefix-defined group."""
    name = header.strip().upper()

    if segment == "VP4":
        # Longest prefixes first.
        for prefix in ("LGWT_", "LGR_", "DG_"):
            if name.startswith(prefix):
                return prefix[:-1]
        return None

    for prefix in ("WTP_", "RP_", "R_", "O_", "G_"):
        if name.startswith(prefix):
            return prefix[:-1]

    return None


def read_fasta(path: Path) -> dict[str, str]:
    """Read an aligned FASTA and verify unique headers/equal alignment length."""
    sequences: dict[str, str] = {}
    header = None
    pieces: list[str] = []

    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    if header in sequences:
                        raise ValueError(
                            f"Duplicate FASTA header in {path.name}: {header}"
                        )
                    sequences[header] = "".join(pieces).replace(" ", "").upper()

                header = line[1:].strip()
                pieces = []
            else:
                pieces.append(line)

        if header is not None:
            if header in sequences:
                raise ValueError(
                    f"Duplicate FASTA header in {path.name}: {header}"
                )
            sequences[header] = "".join(pieces).replace(" ", "").upper()

    if not sequences:
        raise ValueError(f"No sequences found in {path}")

    lengths = {len(seq) for seq in sequences.values()}
    if len(lengths) != 1:
        raise ValueError(
            f"{path.name} is not a valid multiple-sequence alignment: "
            f"sequence lengths differ ({sorted(lengths)})."
        )

    return sequences


def split_groups(
    sequences: dict[str, str],
    segment: str,
) -> tuple[dict[str, list[str]], list[str]]:
    groups_expected = VP4_GROUPS if segment == "VP4" else NON_VP4_GROUPS
    grouped = {group: [] for group in groups_expected}
    unclassified: list[str] = []

    for header in sequences:
        group = classify_header(header, segment)
        if group in grouped:
            grouped[group].append(header)
        else:
            unclassified.append(header)

    # Alphabetical sorting makes outputs deterministic.
    for group in grouped:
        grouped[group].sort()
    unclassified.sort()

    return grouped, unclassified


# ============================================================================
# COORDINATE / STATE HELPERS
# ============================================================================

def ungapped_length(sequence: str) -> int:
    return sum(char != "-" for char in sequence)


def choose_coordinate_reference(
    sequences: dict[str, str],
    segment: str,
) -> tuple[str, str]:
    """
    Deterministically choose the WT coordinate reference.

    Non-VP4 -> R_
    VP4     -> DG_

    Highest ungapped length wins; alphabetical header breaks ties.
    """
    wt_group = "DG" if segment == "VP4" else "R"

    candidates = [
        (header, seq)
        for header, seq in sequences.items()
        if classify_header(header, segment) == wt_group
    ]

    if not candidates:
        raise ValueError(
            f"{segment}: no {wt_group}_ sequence available for WT "
            "coordinate numbering."
        )

    candidates.sort(
        key=lambda item: (-ungapped_length(item[1]), item[0])
    )

    return candidates[0]


def build_position_map(reference_sequence: str) -> list[int | None]:
    """Map alignment columns to 1-based ungapped WT-reference positions."""
    mapping: list[int | None] = []
    position = 0

    for char in reference_sequence:
        if char == "-":
            mapping.append(None)
        else:
            position += 1
            mapping.append(position)

    return mapping


def is_valid_state(char: str, seq_type: str) -> bool:
    if seq_type == "AA":
        return char in AA_VALID
    return char in NU_VALID


def modal_state(chars: list[str], seq_type: str) -> str | None:
    """
    Return the unique modal callable state.

    Tied modes are considered unresolved and return None.
    """
    valid = [char for char in chars if is_valid_state(char, seq_type)]

    if not valid:
        return None

    ranked = Counter(valid).most_common()

    if len(ranked) > 1 and ranked[0][1] == ranked[1][1]:
        return None

    return ranked[0][0]


def group_chars(
    headers: list[str],
    sequences: dict[str, str],
    column: int,
) -> list[str]:
    return [sequences[header][column] for header in headers]


def state_stats(
    headers: list[str],
    sequences: dict[str, str],
    column: int,
    state: str,
    seq_type: str,
) -> dict[str, float | int]:
    """
    Calculate site-specific state frequency using CALLABLE sequences only.

    total_N:
        all sequences belonging to the group

    callable_N:
        sequences with an unambiguous residue/base at this site

    missing_N:
        total_N - callable_N

    n:
        callable sequences carrying the candidate reassortant state

    pct:
        n / callable_N * 100

    callable_fraction:
        callable_N / total_N

    This is the key missing-data correction relative to the earlier script.
    """
    values = group_chars(headers, sequences, column)

    total_N = len(values)
    callable_values = [
        value for value in values
        if is_valid_state(value, seq_type)
    ]

    callable_N = len(callable_values)
    missing_N = total_N - callable_N
    n_state = sum(value == state for value in callable_values)

    pct_value = (
        math.nan
        if callable_N == 0
        else 100.0 * n_state / callable_N
    )

    callable_fraction = (
        math.nan
        if total_N == 0
        else callable_N / total_N
    )

    return {
        "n": n_state,
        "total_N": total_N,
        "callable_N": callable_N,
        "missing_N": missing_N,
        "pct": pct_value,
        "callable_fraction": callable_fraction,
    }


def coverage_is_sufficient(
    stat: dict[str, float | int],
    min_callable_fraction: float,
) -> bool:
    """
    Required existing group has enough callable observations.

    An absent group is handled separately and is not passed here as evidence.
    """
    if stat["total_N"] == 0:
        return False

    if stat["callable_N"] == 0:
        return False

    return float(stat["callable_fraction"]) >= min_callable_fraction


def exclusion_group_passes(
    stat: dict[str, float | int],
    min_callable_fraction: float,
) -> tuple[bool, str]:
    """
    Evaluate an R/O/G or DG exclusion group.

    - Group absent entirely: NA, passes by design.
    - Group exists but too much missing data: unresolved, fails.
    - Adequately callable existing group: reassortant state must be 0%.
    """
    if stat["total_N"] == 0:
        return True, "group_absent_NA"

    if not coverage_is_sufficient(stat, min_callable_fraction):
        return False, "insufficient_callable_coverage"

    if stat["pct"] == 0.0:
        return True, "zero_percent"

    return False, "reassortant_state_present"


# ============================================================================
# OUTPUT SCHEMAS
# ============================================================================

GROUP_ROLES = ["Reassortant", "WT_B", "WT_A", "C", "D"]

SIGNATURE_FIELDS = [
    "Sequence_type",
    "Segment",
    "Position",
    "Alignment_column",
    "Mutation",
    "WT_state",
    "Reassortant_state",
    "Criterion",
    "Coordinate_reference",
    "Min_callable_fraction",
]

for role in GROUP_ROLES:
    SIGNATURE_FIELDS.extend([
        f"{role}_group",
        f"{role}_n",
        f"{role}_total_N",
        f"{role}_callable_N",
        f"{role}_missing_N",
        f"{role}_callable_fraction",
        f"{role}_pct",
        f"{role}_n_over_callable",
    ])

CANDIDATE_AUDIT_FIELDS = SIGNATURE_FIELDS + [
    "Status",
    "Decision_reason",
]

FAILURE_FIELDS = [
    "Sequence_type",
    "Segment",
    "Position",
    "Alignment_column",
    "Candidate_mutation",
    "WT_state",
    "Reassortant_state",
    "Offending_group",
    "Offending_sequence",
    "Observed_state",
    "Expected_state",
    "Failure_type",
    "Decision_reason",
]


# ============================================================================
# ROW BUILDING
# ============================================================================

def blank_stat_fields(role: str) -> dict[str, object]:
    return {
        f"{role}_group": "",
        f"{role}_n": "",
        f"{role}_total_N": "",
        f"{role}_callable_N": "",
        f"{role}_missing_N": "",
        f"{role}_callable_fraction": "",
        f"{role}_pct": "",
        f"{role}_n_over_callable": "",
    }


def add_group_stat(
    row: dict[str, object],
    role: str,
    group_name: str,
    stat: dict[str, float | int],
) -> None:
    row[f"{role}_group"] = group_name
    row[f"{role}_n"] = stat["n"]
    row[f"{role}_total_N"] = stat["total_N"]
    row[f"{role}_callable_N"] = stat["callable_N"]
    row[f"{role}_missing_N"] = stat["missing_N"]

    if stat["total_N"] == 0:
        row[f"{role}_callable_fraction"] = ""
        row[f"{role}_pct"] = ""
        row[f"{role}_n_over_callable"] = "NA"
        return

    row[f"{role}_callable_fraction"] = stat["callable_fraction"]
    row[f"{role}_pct"] = stat["pct"]

    if stat["callable_N"] == 0:
        row[f"{role}_n_over_callable"] = f"{stat['n']}/0"
    else:
        row[f"{role}_n_over_callable"] = (
            f"{stat['n']}/{stat['callable_N']}"
        )


def build_site_row(
    *,
    seq_type: str,
    segment: str,
    position: int,
    alignment_column: int,
    mutation: str,
    wt_state: str,
    reassortant_state: str,
    criterion: str,
    coordinate_reference: str,
    min_callable_fraction: float,
    stats: dict[str, dict[str, float | int]],
) -> dict[str, object]:

    row: dict[str, object] = {
        "Sequence_type": seq_type,
        "Segment": segment,
        "Position": position,
        "Alignment_column": alignment_column,
        "Mutation": mutation,
        "WT_state": wt_state,
        "Reassortant_state": reassortant_state,
        "Criterion": criterion,
        "Coordinate_reference": coordinate_reference,
        "Min_callable_fraction": min_callable_fraction,
    }

    for role in GROUP_ROLES:
        row.update(blank_stat_fields(role))

    if segment == "VP4":
        mapping = {
            "Reassortant": "LGR",
            "WT_B": "LGWT",
            "WT_A": "DG",
        }
    else:
        mapping = {
            "Reassortant": "RP",
            "WT_B": "WTP",
            "WT_A": "R",
            "C": "O",
            "D": "G",
        }

    for role, group_name in mapping.items():
        add_group_stat(row, role, group_name, stats[group_name])

    return row


# ============================================================================
# CORE ANALYSIS
# ============================================================================

def analyze_alignment(
    segment: str,
    seq_type: str,
    fasta_path: Path,
    min_callable_fraction: float,
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    dict[str, list[str]],
    list[str],
]:

    sequences = read_fasta(fasta_path)
    grouped, unclassified = split_groups(sequences, segment)

    groups_expected = VP4_GROUPS if segment == "VP4" else NON_VP4_GROUPS

    if segment == "VP4":
        reassortant_group = "LGR"
        wt_b_group = "LGWT"
        wt_reference_group = "DG"
        exclusion_groups = ["DG"]
    else:
        reassortant_group = "RP"
        wt_b_group = "WTP"
        wt_reference_group = "R"
        exclusion_groups = ["R", "O", "G"]

    print("\n" + "=" * 72)
    print(f"{segment} {seq_type}: {fasta_path.name}")
    print("=" * 72)

    for group in groups_expected:
        print(f"{group:5s}: {len(grouped[group])}")

    if unclassified:
        print(f"UNCLASSIFIED: {len(unclassified)}")

    # Required groups must exist at alignment level.
    missing_required_groups = [
        group
        for group in (reassortant_group, wt_b_group, wt_reference_group)
        if len(grouped[group]) == 0
    ]

    if missing_required_groups:
        print(
            "SKIPPED alignment: missing required group(s): "
            + ", ".join(missing_required_groups)
        )
        return [], [], [], grouped, unclassified

    reference_header, reference_sequence = choose_coordinate_reference(
        sequences,
        segment,
    )
    position_map = build_position_map(reference_sequence)

    print(f"WT coordinate reference: {reference_header}")
    print(
        f"Minimum callable fraction: "
        f"{min_callable_fraction:.2%}"
    )

    accepted: list[dict[str, object]] = []
    candidate_audit: list[dict[str, object]] = []
    failure_details: list[dict[str, object]] = []

    for col in range(len(reference_sequence)):

        position = position_map[col]

        # No WT-reference coordinate at this alignment column.
        if position is None:
            continue

        # --------------------------------------------------------------------
        # Candidate reassortant state
        # --------------------------------------------------------------------

        reassortant_chars = group_chars(
            grouped[reassortant_group],
            sequences,
            col,
        )
        reassortant_state = modal_state(reassortant_chars, seq_type)

        if reassortant_state is None:
            continue

        reass_stat = state_stats(
            grouped[reassortant_group],
            sequences,
            col,
            reassortant_state,
            seq_type,
        )

        # Candidate discovery requires adequate reassortant coverage.
        if not coverage_is_sufficient(
            reass_stat,
            min_callable_fraction,
        ):
            continue

        # Candidate discovery requires 100% among CALLABLE reassortants.
        if reass_stat["pct"] != 100.0:
            continue

        # --------------------------------------------------------------------
        # WT state
        # --------------------------------------------------------------------

        wt_chars = group_chars(
            grouped[wt_reference_group],
            sequences,
            col,
        )
        wt_state = modal_state(wt_chars, seq_type)

        if wt_state is None:
            continue

        if wt_state == reassortant_state:
            continue

        mutation = f"{wt_state}{position}{reassortant_state}"

        # --------------------------------------------------------------------
        # Calculate all group statistics
        # --------------------------------------------------------------------

        stats = {
            group: state_stats(
                grouped[group],
                sequences,
                col,
                reassortant_state,
                seq_type,
            )
            for group in groups_expected
        }

        wt_b_stat = stats[wt_b_group]

        # --------------------------------------------------------------------
        # Evaluate coverage and exclusions
        # --------------------------------------------------------------------

        decision_reasons: list[str] = []

        wt_b_coverage_ok = coverage_is_sufficient(
            wt_b_stat,
            min_callable_fraction,
        )

        if not wt_b_coverage_ok:
            decision_reasons.append(
                f"{wt_b_group} callable coverage "
                f"{int(wt_b_stat['callable_N'])}/"
                f"{int(wt_b_stat['total_N'])} "
                f"({float(wt_b_stat['callable_fraction']):.2%}) "
                f"is below threshold "
                f"{min_callable_fraction:.2%}"
            )

        exclusions_pass = True

        for group in exclusion_groups:
            passed, reason = exclusion_group_passes(
                stats[group],
                min_callable_fraction,
            )

            if not passed:
                exclusions_pass = False

                if reason == "insufficient_callable_coverage":
                    stat = stats[group]
                    decision_reasons.append(
                        f"{group} callable coverage "
                        f"{int(stat['callable_N'])}/"
                        f"{int(stat['total_N'])} "
                        f"({float(stat['callable_fraction']):.2%}) "
                        f"is below threshold "
                        f"{min_callable_fraction:.2%}"
                    )

                elif reason == "reassortant_state_present":
                    stat = stats[group]
                    decision_reasons.append(
                        f"{group} contains reassortant state "
                        f"{reassortant_state} in "
                        f"{int(stat['n'])}/"
                        f"{int(stat['callable_N'])} callable sequences "
                        f"({float(stat['pct']):.2f}%)"
                    )

        # --------------------------------------------------------------------
        # Apply the two predefined criteria
        # --------------------------------------------------------------------

        criterion1 = (
            wt_b_coverage_ok
            and wt_b_stat["pct"] == 100.0
            and exclusions_pass
        )

        criterion2 = (
            wt_b_coverage_ok
            and not math.isnan(float(wt_b_stat["pct"]))
            and 0.0 < float(wt_b_stat["pct"]) < 100.0
            and exclusions_pass
        )

        if criterion1:
            criterion = "Criterion_1_Strict_B_signature"
            status = "ACCEPTED"
            decision_reason = (
                "Reassortant group fixed at 100% among callable sequences; "
                f"{wt_b_group}=100%; all existing exclusion groups=0%; "
                "all required existing groups met callable-data threshold."
            )

        elif criterion2:
            criterion = "Criterion_2_Reassortant_enriched_B_signature"
            status = "ACCEPTED"
            decision_reason = (
                "Reassortant group fixed at 100% among callable sequences; "
                f"0%<{wt_b_group}<100%; all existing exclusion groups=0%; "
                "all required existing groups met callable-data threshold."
            )

        else:
            criterion = ""
            status = "REJECTED"

            if (
                wt_b_coverage_ok
                and wt_b_stat["pct"] == 0.0
            ):
                decision_reasons.append(
                    f"{wt_b_group}=0%, so site satisfies neither "
                    "Criterion 1 nor Criterion 2."
                )

            if not decision_reasons:
                decision_reasons.append(
                    "Candidate did not satisfy either complete criterion."
                )

            decision_reason = "; ".join(decision_reasons)

        site_row = build_site_row(
            seq_type=seq_type,
            segment=segment,
            position=position,
            alignment_column=col + 1,
            mutation=mutation,
            wt_state=wt_state,
            reassortant_state=reassortant_state,
            criterion=criterion,
            coordinate_reference=reference_header,
            min_callable_fraction=min_callable_fraction,
            stats=stats,
        )

        audit_row = dict(site_row)
        audit_row["Status"] = status
        audit_row["Decision_reason"] = decision_reason
        candidate_audit.append(audit_row)

        if status == "ACCEPTED":
            accepted.append(site_row)
        else:
            failure_details.extend(
                build_failure_details(
                    seq_type=seq_type,
                    segment=segment,
                    position=position,
                    alignment_column=col + 1,
                    mutation=mutation,
                    wt_state=wt_state,
                    reassortant_state=reassortant_state,
                    decision_reason=decision_reason,
                    grouped=grouped,
                    sequences=sequences,
                    column=col,
                    seq_type_for_calls=seq_type,
                    wt_b_group=wt_b_group,
                    exclusion_groups=exclusion_groups,
                    stats=stats,
                    min_callable_fraction=min_callable_fraction,
                )
            )

    print(f"Accepted signatures: {len(accepted)}")
    print(f"Candidate sites audited: {len(candidate_audit)}")

    return (
        accepted,
        candidate_audit,
        failure_details,
        grouped,
        unclassified,
    )


def build_failure_details(
    *,
    seq_type: str,
    segment: str,
    position: int,
    alignment_column: int,
    mutation: str,
    wt_state: str,
    reassortant_state: str,
    decision_reason: str,
    grouped: dict[str, list[str]],
    sequences: dict[str, str],
    column: int,
    seq_type_for_calls: str,
    wt_b_group: str,
    exclusion_groups: list[str],
    stats: dict[str, dict[str, float | int]],
    min_callable_fraction: float,
) -> list[dict[str, object]]:

    rows: list[dict[str, object]] = []

    # WTP/LGWT: list callable sequences that differ from reassortant state.
    for header in grouped[wt_b_group]:
        observed = sequences[header][column]

        if not is_valid_state(observed, seq_type_for_calls):
            # Missing calls are documented explicitly.
            rows.append({
                "Sequence_type": seq_type,
                "Segment": segment,
                "Position": position,
                "Alignment_column": alignment_column,
                "Candidate_mutation": mutation,
                "WT_state": wt_state,
                "Reassortant_state": reassortant_state,
                "Offending_group": wt_b_group,
                "Offending_sequence": header,
                "Observed_state": observed,
                "Expected_state": reassortant_state,
                "Failure_type": "Missing/non-callable state",
                "Decision_reason": decision_reason,
            })

        elif observed != reassortant_state:
            rows.append({
                "Sequence_type": seq_type,
                "Segment": segment,
                "Position": position,
                "Alignment_column": alignment_column,
                "Candidate_mutation": mutation,
                "WT_state": wt_state,
                "Reassortant_state": reassortant_state,
                "Offending_group": wt_b_group,
                "Offending_sequence": header,
                "Observed_state": observed,
                "Expected_state": reassortant_state,
                "Failure_type": "Callable state differs from reassortant state",
                "Decision_reason": decision_reason,
            })

    # Exclusion groups: list missing calls and any callable reassortant-state calls.
    for group in exclusion_groups:
        if not grouped[group]:
            continue

        stat = stats[group]
        coverage_low = not coverage_is_sufficient(
            stat,
            min_callable_fraction,
        )

        for header in grouped[group]:
            observed = sequences[header][column]

            if not is_valid_state(observed, seq_type_for_calls):
                if coverage_low:
                    rows.append({
                        "Sequence_type": seq_type,
                        "Segment": segment,
                        "Position": position,
                        "Alignment_column": alignment_column,
                        "Candidate_mutation": mutation,
                        "WT_state": wt_state,
                        "Reassortant_state": reassortant_state,
                        "Offending_group": group,
                        "Offending_sequence": header,
                        "Observed_state": observed,
                        "Expected_state": f"callable non-{reassortant_state}",
                        "Failure_type": "Missing/non-callable state contributes to low coverage",
                        "Decision_reason": decision_reason,
                    })

            elif observed == reassortant_state:
                rows.append({
                    "Sequence_type": seq_type,
                    "Segment": segment,
                    "Position": position,
                    "Alignment_column": alignment_column,
                    "Candidate_mutation": mutation,
                    "WT_state": wt_state,
                    "Reassortant_state": reassortant_state,
                    "Offending_group": group,
                    "Offending_sequence": header,
                    "Observed_state": observed,
                    "Expected_state": f"not {reassortant_state}",
                    "Failure_type": "Exclusion group carries reassortant state",
                    "Decision_reason": decision_reason,
                })

    if not rows:
        rows.append({
            "Sequence_type": seq_type,
            "Segment": segment,
            "Position": position,
            "Alignment_column": alignment_column,
            "Candidate_mutation": mutation,
            "WT_state": wt_state,
            "Reassortant_state": reassortant_state,
            "Offending_group": "",
            "Offending_sequence": "",
            "Observed_state": "",
            "Expected_state": "",
            "Failure_type": "Site-level criterion failure",
            "Decision_reason": decision_reason,
        })

    return rows


# ============================================================================
# CSV WRITING
# ============================================================================

def format_output_value(value: object) -> object:
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.6f}"
    return value


def write_dict_rows(
    path: Path,
    rows: list[dict[str, object]],
    fields: list[str],
) -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()

        for row in rows:
            formatted = {
                key: format_output_value(row.get(key, ""))
                for key in fields
            }
            writer.writerow(formatted)


def write_group_counts(
    path: Path,
    rows: list[dict[str, object]],
) -> None:
    fields = [
        "Segment",
        "Sequence_type",
        "Group",
        "N_sequences",
        "Unclassified_N",
    ]
    write_dict_rows(path, rows, fields)


# ============================================================================
# GRAPHING
# ============================================================================

def numeric(value: object) -> float:
    if value in ("", None):
        return np.nan
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def make_signature_heatmap(
    rows: list[dict[str, object]],
    seq_type: str,
    output_dir: Path,
) -> None:
    """Heatmap of reassortant-state prevalence among callable sequences."""
    subset = [
        row for row in rows
        if row["Sequence_type"] == seq_type
    ]

    if not subset:
        print(f"No accepted {seq_type} signatures; prevalence heatmap skipped.")
        return

    segment_order = {
        segment: i for i, segment in enumerate(SEGMENTS)
    }

    subset.sort(
        key=lambda row: (
            segment_order[str(row["Segment"])],
            int(row["Position"]),
            str(row["Mutation"]),
        )
    )

    labels = [
        f"{row['Segment']} {row['Mutation']}"
        for row in subset
    ]

    roles = ["Reassortant", "WT_B", "WT_A", "C", "D"]
    col_labels = [
        "Reassortant\nRP/LGR",
        "WT-B/light-green\nWTP/LGWT",
        "WT-A/deep-green\nR/DG",
        "Type C\nO",
        "Type D\nG",
    ]

    matrix = np.array(
        [
            [numeric(row[f"{role}_pct"]) for role in roles]
            for row in subset
        ],
        dtype=float,
    )

    fig_height = max(4.5, 0.38 * len(labels) + 2.5)
    fig, ax = plt.subplots(figsize=(10, fig_height))

    masked = np.ma.masked_invalid(matrix)
    image = ax.imshow(
        masked,
        aspect="auto",
        vmin=0,
        vmax=100,
    )

    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels)
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_title(
        f"{seq_type} reassortant-state prevalence among callable sequences"
    )

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            text = "NA" if np.isnan(value) else f"{value:.0f}%"
            ax.text(
                j, i, text,
                ha="center",
                va="center",
                fontsize=8,
            )

    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("Reassortant-state prevalence (%)")

    plt.tight_layout()
    plt.savefig(
        output_dir / f"P1_{seq_type}_signature_prevalence_heatmap.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.savefig(
        output_dir / f"P1_{seq_type}_signature_prevalence_heatmap.pdf",
        bbox_inches="tight",
    )
    plt.close()


def make_coverage_heatmap(
    rows: list[dict[str, object]],
    seq_type: str,
    output_dir: Path,
) -> None:
    """Heatmap showing callable-data fraction for accepted signatures."""
    subset = [
        row for row in rows
        if row["Sequence_type"] == seq_type
    ]

    if not subset:
        return

    segment_order = {
        segment: i for i, segment in enumerate(SEGMENTS)
    }

    subset.sort(
        key=lambda row: (
            segment_order[str(row["Segment"])],
            int(row["Position"]),
            str(row["Mutation"]),
        )
    )

    labels = [
        f"{row['Segment']} {row['Mutation']}"
        for row in subset
    ]

    roles = ["Reassortant", "WT_B", "WT_A", "C", "D"]
    col_labels = [
        "Reassortant\nRP/LGR",
        "WT-B/light-green\nWTP/LGWT",
        "WT-A/deep-green\nR/DG",
        "Type C\nO",
        "Type D\nG",
    ]

    matrix = np.array(
        [
            [
                100.0 * numeric(row[f"{role}_callable_fraction"])
                for role in roles
            ]
            for row in subset
        ],
        dtype=float,
    )

    fig_height = max(4.5, 0.38 * len(labels) + 2.5)
    fig, ax = plt.subplots(figsize=(10, fig_height))

    masked = np.ma.masked_invalid(matrix)
    image = ax.imshow(
        masked,
        aspect="auto",
        vmin=0,
        vmax=100,
    )

    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels)
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_title(
        f"{seq_type} callable-data coverage at accepted signature sites"
    )

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            text = "NA" if np.isnan(value) else f"{value:.0f}%"
            ax.text(
                j, i, text,
                ha="center",
                va="center",
                fontsize=8,
            )

    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("Callable sequences (%)")

    plt.tight_layout()
    plt.savefig(
        output_dir / f"P1_{seq_type}_signature_coverage_heatmap.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.savefig(
        output_dir / f"P1_{seq_type}_signature_coverage_heatmap.pdf",
        bbox_inches="tight",
    )
    plt.close()



# ============================================================================
# FISHER'S EXACT TESTS FOR ACCEPTED SIGNATURE SITES
# ============================================================================

FISHER_FIELDS = [
    "Sequence_type",
    "Segment",
    "Position",
    "Alignment_column",
    "Mutation",
    "Criterion",
    "Reassortant_state",
    "Test_mode",
    "Reassortant_groups",
    "Comparison_group",
    "Comparison_role",
    "Reassortant_state_n",
    "Reassortant_callable_N",
    "Reassortant_not_state_n",
    "Comparison_state_n",
    "Comparison_callable_N",
    "Comparison_not_state_n",
    "Reassortant_n_over_callable",
    "Comparison_n_over_callable",
    "Odds_ratio",
    "Fisher_p_two_sided",
    "BH_FDR_within_family",
    "BH_FDR_global",
    "Significant_FDR_0.05",
]


def bh_fdr(p_values: list[float]) -> list[float]:
    """
    Benjamini-Hochberg false-discovery-rate correction.

    NaN p-values are retained as NaN and are not included in the correction.
    """
    adjusted = [math.nan] * len(p_values)

    valid = [
        (index, float(p))
        for index, p in enumerate(p_values)
        if p is not None and not math.isnan(float(p))
    ]

    if not valid:
        return adjusted

    valid.sort(key=lambda item: item[1])
    m = len(valid)

    ranked_adjusted = [0.0] * m

    for rank, (_, p_value) in enumerate(valid, start=1):
        ranked_adjusted[rank - 1] = min(
            1.0,
            p_value * m / rank,
        )

    # Enforce monotonicity from largest rank to smallest.
    for i in range(m - 2, -1, -1):
        ranked_adjusted[i] = min(
            ranked_adjusted[i],
            ranked_adjusted[i + 1],
        )

    for (original_index, _), adj_p in zip(valid, ranked_adjusted):
        adjusted[original_index] = adj_p

    return adjusted


def _int_field(row: dict[str, object], field: str) -> int:
    value = row.get(field, "")

    if value in ("", None):
        return 0

    return int(float(value))


def _fisher_row(
    *,
    signature_row: dict[str, object],
    test_mode: str,
    reassortant_groups: str,
    comparison_group: str,
    comparison_role: str,
    reassortant_state_n: int,
    reassortant_callable_n: int,
    comparison_state_n: int,
    comparison_callable_n: int,
) -> dict[str, object] | None:
    """
    Build and run one 2x2 Fisher's exact test.

    Table:

                         reassortant state   other callable state
        test group              a                    b
        comparison group        c                    d
    """
    if reassortant_callable_n <= 0 or comparison_callable_n <= 0:
        return None

    a = reassortant_state_n
    b = reassortant_callable_n - reassortant_state_n
    c = comparison_state_n
    d = comparison_callable_n - comparison_state_n

    odds_ratio, p_value = fisher_exact(
        [[a, b], [c, d]],
        alternative="two-sided",
    )

    return {
        "Sequence_type": signature_row["Sequence_type"],
        "Segment": signature_row["Segment"],
        "Position": signature_row["Position"],
        "Alignment_column": signature_row["Alignment_column"],
        "Mutation": signature_row["Mutation"],
        "Criterion": signature_row["Criterion"],
        "Reassortant_state": signature_row["Reassortant_state"],
        "Test_mode": test_mode,
        "Reassortant_groups": reassortant_groups,
        "Comparison_group": comparison_group,
        "Comparison_role": comparison_role,
        "Reassortant_state_n": a,
        "Reassortant_callable_N": reassortant_callable_n,
        "Reassortant_not_state_n": b,
        "Comparison_state_n": c,
        "Comparison_callable_N": comparison_callable_n,
        "Comparison_not_state_n": d,
        "Reassortant_n_over_callable": f"{a}/{reassortant_callable_n}",
        "Comparison_n_over_callable": f"{c}/{comparison_callable_n}",
        "Odds_ratio": float(odds_ratio),
        "Fisher_p_two_sided": float(p_value),
        "BH_FDR_within_family": math.nan,
        "BH_FDR_global": math.nan,
        "Significant_FDR_0.05": "",
    }


def run_fisher_tests(
    signatures: list[dict[str, object]],
) -> list[dict[str, object]]:
    """
    Run Fisher's exact tests at every ACCEPTED Criterion 1 or Criterion 2 site.

    Two predefined test modes are evaluated.

    Non-VP4:
        1. RP only versus each existing comparison group R, O, G.
        2. RP + WTP combined versus each existing comparison group R, O, G.

    VP4:
        1. LGR only versus DG.
        2. LGR + LGWT combined versus DG.

    Only callable sequences enter the 2x2 tables. Missing/gapped/ambiguous
    calls remain excluded exactly as in the signature-discovery analysis.
    """
    fisher_rows: list[dict[str, object]] = []

    for row in signatures:
        segment = str(row["Segment"])

        reass_n = _int_field(row, "Reassortant_n")
        reass_callable = _int_field(row, "Reassortant_callable_N")

        wt_b_n = _int_field(row, "WT_B_n")
        wt_b_callable = _int_field(row, "WT_B_callable_N")

        if segment == "VP4":
            comparison_specs = [
                ("WT_A", str(row.get("WT_A_group", "DG"))),
            ]

            single_label = "LGR_only"
            single_groups = "LGR"
            combined_label = "LGR_plus_LGWT"
            combined_groups = "LGR+LGWT"

        else:
            comparison_specs = [
                ("WT_A", str(row.get("WT_A_group", "R"))),
                ("C", str(row.get("C_group", "O"))),
                ("D", str(row.get("D_group", "G"))),
            ]

            single_label = "RP_only"
            single_groups = "RP"
            combined_label = "RP_plus_WTP"
            combined_groups = "RP+WTP"

        for comparison_role, comparison_group in comparison_specs:
            comparison_n = _int_field(
                row,
                f"{comparison_role}_n",
            )
            comparison_callable = _int_field(
                row,
                f"{comparison_role}_callable_N",
            )

            # Group absent at alignment level -> no Fisher test.
            if comparison_callable <= 0:
                continue

            # ----------------------------------------------------------
            # Test 1: reassortant group only vs comparison group
            # ----------------------------------------------------------
            single = _fisher_row(
                signature_row=row,
                test_mode=single_label,
                reassortant_groups=single_groups,
                comparison_group=comparison_group,
                comparison_role=comparison_role,
                reassortant_state_n=reass_n,
                reassortant_callable_n=reass_callable,
                comparison_state_n=comparison_n,
                comparison_callable_n=comparison_callable,
            )

            if single is not None:
                fisher_rows.append(single)

            # ----------------------------------------------------------
            # Test 2: reassortant + WT-B group combined vs comparison
            # ----------------------------------------------------------
            combined = _fisher_row(
                signature_row=row,
                test_mode=combined_label,
                reassortant_groups=combined_groups,
                comparison_group=comparison_group,
                comparison_role=comparison_role,
                reassortant_state_n=reass_n + wt_b_n,
                reassortant_callable_n=reass_callable + wt_b_callable,
                comparison_state_n=comparison_n,
                comparison_callable_n=comparison_callable,
            )

            if combined is not None:
                fisher_rows.append(combined)

    # --------------------------------------------------------------------
    # Benjamini-Hochberg correction WITHIN each predefined test family.
    #
    # A family is:
    #   sequence type + test mode + comparison group
    #
    # Example:
    #   AA + RP_only + R
    #   AA + RP_plus_WTP + R
    #   NU + RP_only + O
    #
    # This corrects across all accepted sites tested for the same biological
    # comparison while keeping distinct comparisons separate.
    # --------------------------------------------------------------------
    families: dict[tuple[str, str, str], list[int]] = {}

    for index, row in enumerate(fisher_rows):
        key = (
            str(row["Sequence_type"]),
            str(row["Test_mode"]),
            str(row["Comparison_group"]),
        )
        families.setdefault(key, []).append(index)

    for indices in families.values():
        p_values = [
            float(fisher_rows[index]["Fisher_p_two_sided"])
            for index in indices
        ]

        adjusted = bh_fdr(p_values)

        for index, adj_p in zip(indices, adjusted):
            fisher_rows[index]["BH_FDR_within_family"] = adj_p

    # Also report a conservative BH correction across every Fisher test.
    all_p_values = [
        float(row["Fisher_p_two_sided"])
        for row in fisher_rows
    ]
    global_adjusted = bh_fdr(all_p_values)

    for row, global_adj in zip(fisher_rows, global_adjusted):
        row["BH_FDR_global"] = global_adj

        family_adj = float(row["BH_FDR_within_family"])

        row["Significant_FDR_0.05"] = (
            "YES"
            if not math.isnan(family_adj) and family_adj < 0.05
            else "NO"
        )

    segment_order = {
        segment: i
        for i, segment in enumerate(SEGMENTS)
    }
    type_order = {"AA": 0, "NU": 1}

    fisher_rows.sort(
        key=lambda row: (
            type_order[str(row["Sequence_type"])],
            segment_order[str(row["Segment"])],
            int(row["Position"]),
            str(row["Test_mode"]),
            str(row["Comparison_group"]),
        )
    )

    return fisher_rows


def make_fisher_fdr_heatmap(
    rows: list[dict[str, object]],
    seq_type: str,
    output_dir: Path,
) -> None:
    """
    Heatmap of -log10(BH-FDR) for Fisher tests at accepted signature sites.

    The plotted FDR is the within-family Benjamini-Hochberg adjusted p-value.
    """
    subset = [
        row for row in rows
        if row["Sequence_type"] == seq_type
    ]

    if not subset:
        print(f"No {seq_type} Fisher-test rows; Fisher heatmap skipped.")
        return

    segment_order = {
        segment: i
        for i, segment in enumerate(SEGMENTS)
    }

    # Unique signature sites.
    site_keys = sorted(
        {
            (
                str(row["Segment"]),
                int(row["Position"]),
                str(row["Mutation"]),
                str(row["Criterion"]),
            )
            for row in subset
        },
        key=lambda item: (
            segment_order[item[0]],
            item[1],
            item[2],
        ),
    )

    # Unique test columns, deterministic and biologically readable.
    preferred_test_order = [
        "RP_only vs R",
        "RP_only vs O",
        "RP_only vs G",
        "RP_plus_WTP vs R",
        "RP_plus_WTP vs O",
        "RP_plus_WTP vs G",
        "LGR_only vs DG",
        "LGR_plus_LGWT vs DG",
    ]

    observed_tests = {
        f"{row['Test_mode']} vs {row['Comparison_group']}"
        for row in subset
    }

    test_labels = [
        label
        for label in preferred_test_order
        if label in observed_tests
    ]

    # Include any unexpected but valid labels deterministically.
    test_labels.extend(
        sorted(observed_tests - set(test_labels))
    )

    site_index = {
        key: i
        for i, key in enumerate(site_keys)
    }
    test_index = {
        label: j
        for j, label in enumerate(test_labels)
    }

    matrix = np.full(
        (len(site_keys), len(test_labels)),
        np.nan,
        dtype=float,
    )

    for row in subset:
        key = (
            str(row["Segment"]),
            int(row["Position"]),
            str(row["Mutation"]),
            str(row["Criterion"]),
        )

        label = (
            f"{row['Test_mode']} vs "
            f"{row['Comparison_group']}"
        )

        fdr = float(row["BH_FDR_within_family"])

        if math.isnan(fdr):
            value = np.nan
        elif fdr <= 0:
            # Avoid infinite plotting values if a numerical p-value is exactly zero.
            value = 300.0
        else:
            value = -math.log10(fdr)

        matrix[
            site_index[key],
            test_index[label],
        ] = value

    labels = [
        f"{segment} {mutation}"
        for segment, _, mutation, _ in site_keys
    ]

    fig_height = max(
        5.0,
        0.34 * len(labels) + 2.5,
    )
    fig_width = max(
        10.0,
        1.5 * len(test_labels) + 4.0,
    )

    fig, ax = plt.subplots(
        figsize=(fig_width, fig_height)
    )

    masked = np.ma.masked_invalid(matrix)

    image = ax.imshow(
        masked,
        aspect="auto",
    )

    ax.set_xticks(
        np.arange(len(test_labels))
    )
    ax.set_xticklabels(
        test_labels,
        rotation=45,
        ha="right",
    )

    ax.set_yticks(
        np.arange(len(labels))
    )
    ax.set_yticklabels(labels)

    ax.set_title(
        f"{seq_type} signature-site Fisher tests "
        "(-log10 BH-FDR)"
    )

    threshold = -math.log10(0.05)

    # Annotate each populated cell with adjusted p-value.
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if np.isnan(matrix[i, j]):
                text = "NA"
            else:
                # Find original FDR for display.
                site_key = site_keys[i]
                test_label = test_labels[j]

                matched = [
                    row for row in subset
                    if (
                        str(row["Segment"]),
                        int(row["Position"]),
                        str(row["Mutation"]),
                        str(row["Criterion"]),
                    ) == site_key
                    and (
                        f"{row['Test_mode']} vs "
                        f"{row['Comparison_group']}"
                    ) == test_label
                ]

                if matched:
                    fdr = float(
                        matched[0]["BH_FDR_within_family"]
                    )
                    text = (
                        "<1e-300"
                        if fdr == 0
                        else f"{fdr:.1e}"
                    )
                else:
                    text = "NA"

            ax.text(
                j,
                i,
                text,
                ha="center",
                va="center",
                fontsize=7,
            )

    cbar = fig.colorbar(
        image,
        ax=ax,
    )
    cbar.set_label(
        "-log10(BH-FDR)"
    )

    # Explain the significance threshold in a compact caption-like note.
    ax.text(
        1.0,
        1.01,
        f"FDR < 0.05 corresponds to -log10(FDR) > {threshold:.2f}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
    )

    plt.tight_layout()

    plt.savefig(
        output_dir
        / f"P1_{seq_type}_signature_Fisher_BH_FDR_heatmap.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.savefig(
        output_dir
        / f"P1_{seq_type}_signature_Fisher_BH_FDR_heatmap.pdf",
        bbox_inches="tight",
    )

    plt.close()


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    all_signatures: list[dict[str, object]] = []
    all_candidates: list[dict[str, object]] = []
    all_failures: list[dict[str, object]] = []
    group_count_rows: list[dict[str, object]] = []
    input_files_used: list[Path] = []

    for segment in SEGMENTS:
        for seq_type in SEQ_TYPES:

            fasta_path = args.input_dir / f"{segment}_{seq_type}.fasta"

            if not fasta_path.exists():
                print(f"Not found: {fasta_path.name} — skipped.")
                continue

            input_files_used.append(fasta_path)

            (
                signatures,
                candidates,
                failures,
                grouped,
                unclassified,
            ) = analyze_alignment(
                segment=segment,
                seq_type=seq_type,
                fasta_path=fasta_path,
                min_callable_fraction=args.min_callable_fraction,
            )

            all_signatures.extend(signatures)
            all_candidates.extend(candidates)
            all_failures.extend(failures)

            for group, headers in grouped.items():
                group_count_rows.append({
                    "Segment": segment,
                    "Sequence_type": seq_type,
                    "Group": group,
                    "N_sequences": len(headers),
                    "Unclassified_N": len(unclassified),
                })

    # Deterministic sort order.
    segment_order = {
        segment: i for i, segment in enumerate(SEGMENTS)
    }
    type_order = {"AA": 0, "NU": 1}

    def signature_sort_key(row: dict[str, object]):
        return (
            type_order[str(row["Sequence_type"])],
            segment_order[str(row["Segment"])],
            int(row["Position"]),
            str(row["Mutation"]),
        )

    all_signatures.sort(key=signature_sort_key)
    all_candidates.sort(key=signature_sort_key)

    all_failures.sort(
        key=lambda row: (
            type_order[str(row["Sequence_type"])],
            segment_order[str(row["Segment"])],
            int(row["Position"]),
            str(row["Candidate_mutation"]),
            str(row["Offending_group"]),
            str(row["Offending_sequence"]),
        )
    )

    # ------------------------------------------------------------------------
    # Write accepted-signature tables
    # ------------------------------------------------------------------------

    write_dict_rows(
        args.output_dir / "P1_signature_changes.csv",
        all_signatures,
        SIGNATURE_FIELDS,
    )

    write_dict_rows(
        args.output_dir / "P1_AA_signature_changes.csv",
        [
            row for row in all_signatures
            if row["Sequence_type"] == "AA"
        ],
        SIGNATURE_FIELDS,
    )

    write_dict_rows(
        args.output_dir / "P1_NU_signature_changes.csv",
        [
            row for row in all_signatures
            if row["Sequence_type"] == "NU"
        ],
        SIGNATURE_FIELDS,
    )

    # Every candidate fixed in RP/LGR after coverage filtering, accepted or rejected.
    write_dict_rows(
        args.output_dir / "P1_signature_candidate_audit.csv",
        all_candidates,
        CANDIDATE_AUDIT_FIELDS,
    )

    # Sequence-level details explaining rejected candidate sites.
    write_dict_rows(
        args.output_dir / "P1_signature_failure_details.csv",
        all_failures,
        FAILURE_FIELDS,
    )

    write_group_counts(
        args.output_dir / "P1_signature_group_counts.csv",
        group_count_rows,
    )

    write_run_metadata(
        output_dir=args.output_dir,
        input_files=input_files_used,
        min_callable_fraction=args.min_callable_fraction,
    )

    # ------------------------------------------------------------------------
    # Fisher's exact tests for every accepted Criterion 1 / Criterion 2 site
    # ------------------------------------------------------------------------

    fisher_rows = run_fisher_tests(
        all_signatures
    )

    write_dict_rows(
        args.output_dir / "P1_signature_Fisher_tests.csv",
        fisher_rows,
        FISHER_FIELDS,
    )

    write_dict_rows(
        args.output_dir / "P1_AA_signature_Fisher_tests.csv",
        [
            row for row in fisher_rows
            if row["Sequence_type"] == "AA"
        ],
        FISHER_FIELDS,
    )

    write_dict_rows(
        args.output_dir / "P1_NU_signature_Fisher_tests.csv",
        [
            row for row in fisher_rows
            if row["Sequence_type"] == "NU"
        ],
        FISHER_FIELDS,
    )

    for seq_type in SEQ_TYPES:
        make_fisher_fdr_heatmap(
            fisher_rows,
            seq_type,
            args.output_dir,
        )

    # ------------------------------------------------------------------------
    # Figures
    # ------------------------------------------------------------------------

    for seq_type in SEQ_TYPES:
        make_signature_heatmap(
            all_signatures,
            seq_type,
            args.output_dir,
        )
        make_coverage_heatmap(
            all_signatures,
            seq_type,
            args.output_dir,
        )

    # ------------------------------------------------------------------------
    # Console summary
    # ------------------------------------------------------------------------

    aa_n = sum(
        row["Sequence_type"] == "AA"
        for row in all_signatures
    )
    nu_n = sum(
        row["Sequence_type"] == "NU"
        for row in all_signatures
    )

    print("\n" + "=" * 72)
    print("ANALYSIS COMPLETE")
    print("=" * 72)
    print(f"Script version:                {SCRIPT_VERSION}")
    print(f"Minimum callable fraction:    {args.min_callable_fraction:.2%}")
    print(f"Accepted AA signatures:       {aa_n}")
    print(f"Accepted NU signatures:       {nu_n}")
    print(f"Candidate sites audited:      {len(all_candidates)}")
    print(f"Sequence-level failure rows:  {len(all_failures)}")
    print(f"Fisher exact-test rows:       {len(fisher_rows)}")
    print(f"Output directory:             {args.output_dir}")

    print("\nKey missing-data rule:")
    print(
        "  Percentages use callable_N as the denominator. Missing/ambiguous "
        "calls are excluded from the percentage but explicitly reported."
    )
    print(
        "  Existing required groups must meet the minimum callable fraction; "
        "groups absent entirely are reported as NA."
    )


if __name__ == "__main__":
    main()
