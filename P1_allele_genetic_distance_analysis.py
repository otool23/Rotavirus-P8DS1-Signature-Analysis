#!/usr/bin/env python3
"""
Paper 1: reproducible allele genetic-distance analysis
======================================================

Purpose
-------
Quantify within- and between-allele genetic differentiation for Rotavirus A
segments using pairwise percent-identity matrices. The analysis compares allele
B only against allele A, C, or D when those groups are present.

Expected input files
--------------------
VP1.csv, VP2.csv, VP3.csv, VP4.csv, VP6.csv,
NSP1.csv, NSP2.csv, NSP3.csv, NSP4.csv

VP7 and NSP5 are intentionally excluded.

Input matrix format
-------------------
Each CSV must be a square pairwise percent-identity matrix whose first row and
first column contain the same sequence labels in the same order.

Allele labels
-------------
Each sequence label must begin with A-, B-, C-, or D-, for example:
A-KJ123456, B-KC443369, C-LC123456, D-MN123456.

Calculations
------------
p-distance (%) = 100 - percent identity

d_net = d_between - (d_within_B + d_within_other) / 2

Identity ratio = within-B identity / between-group identity
Distance ratio = within-B distance / between-group distance

Interpretation
--------------
Identity ratio > 1 indicates that allele B is more similar internally than it
is to the comparison allele. Distance ratio < 1 indicates that within-B
genetic distance is smaller than B-to-other genetic distance.

Reproducibility
---------------
The script validates matrix shape, label order, symmetry, and diagonal values;
records Python/NumPy/Matplotlib versions and SHA-256 checksums; and writes
summary tables plus PNG/PDF figures. No random sampling is used.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import platform
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_VERSION = "1.0.0"
SEGMENTS = ["VP1", "VP2", "VP3", "VP4", "VP6", "NSP1", "NSP2", "NSP3", "NSP4"]
COMPARISON_ORDER = ["A", "C", "D"]
ALLELE_PATTERN = re.compile(r"^\s*([ABCD])-", flags=re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Reproducible allele genetic-distance analysis.")
    parser.add_argument("--input-dir", type=Path, default=script_dir,
                        help="Directory containing VP1.csv ... NSP4.csv (default: script directory).")
    parser.add_argument("--output-dir", type=Path, default=script_dir / "P1_distance_results",
                        help="Directory for result tables and figures.")
    parser.add_argument("--symmetry-tolerance", type=float, default=1e-6)
    parser.add_argument("--diagonal-tolerance", type=float, default=1e-6)
    args = parser.parse_args()
    args.input_dir = args.input_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    return args


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_metadata(output_dir: Path, input_files: list[Path], args: argparse.Namespace) -> None:
    lines = [
        "Paper 1 allele genetic-distance analysis",
        "========================================",
        f"Script version: {SCRIPT_VERSION}",
        f"Run timestamp (UTC): {datetime.now(timezone.utc).isoformat()}",
        f"Python version: {platform.python_version()}",
        f"Python executable: {sys.executable}",
        f"Platform: {platform.platform()}",
        f"NumPy version: {np.__version__}",
        f"Matplotlib version: {matplotlib.__version__}",
        f"Symmetry tolerance: {args.symmetry_tolerance}",
        f"Diagonal tolerance: {args.diagonal_tolerance}",
        "",
        "Input file SHA-256 checksums",
        "----------------------------",
    ]
    for path in sorted(input_files, key=lambda p: p.name):
        lines.append(f"{path.name}\t{sha256_file(path)}")
    (output_dir / "P1_distance_run_metadata.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def allele(label: str) -> str | None:
    match = ALLELE_PATTERN.match(str(label))
    return match.group(1).upper() if match else None


def parse_numeric(value: str) -> float:
    value = value.strip()
    if value in ("", "NA", "NaN", "nan", "NAN"):
        return math.nan
    return float(value)


def read_identity_matrix(path: Path, symmetry_tolerance: float, diagonal_tolerance: float):
    validation_messages = []
    with open(path, newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.reader(handle))
    if len(rows) < 2 or len(rows[0]) < 2:
        raise ValueError(f"{path.name}: invalid matrix structure.")

    column_labels = [x.strip() for x in rows[0][1:]]
    row_labels = [r[0].strip() for r in rows[1:]]
    matrix = []
    for row_number, row in enumerate(rows[1:], start=2):
        if len(row) != len(column_labels) + 1:
            raise ValueError(f"{path.name}: row {row_number} has wrong number of values.")
        matrix.append([parse_numeric(x) for x in row[1:]])

    if len(row_labels) != len(column_labels):
        raise ValueError(f"{path.name}: matrix is not square.")
    if row_labels != column_labels:
        raise ValueError(f"{path.name}: row and column labels are not identical and in the same order.")

    unclassified = [label for label in row_labels if allele(label) is None]
    if unclassified:
        validation_messages.append(
            f"Unclassified labels: {len(unclassified)}; examples: {', '.join(unclassified[:5])}"
        )

    diagonal_warnings = 0
    for i in range(len(row_labels)):
        v = matrix[i][i]
        if math.isnan(v) or abs(v - 100.0) > diagonal_tolerance:
            diagonal_warnings += 1
    if diagonal_warnings:
        validation_messages.append(
            f"Diagonal warnings: {diagonal_warnings} entries are missing or differ from 100% identity."
        )

    symmetry_warnings = 0
    max_asymmetry = 0.0
    for i in range(len(row_labels)):
        for j in range(i + 1, len(row_labels)):
            a, b = matrix[i][j], matrix[j][i]
            if math.isnan(a) or math.isnan(b):
                continue
            diff = abs(a - b)
            max_asymmetry = max(max_asymmetry, diff)
            if diff > symmetry_tolerance:
                symmetry_warnings += 1
    if symmetry_warnings:
        validation_messages.append(
            f"Symmetry warnings: {symmetry_warnings} pair(s); maximum asymmetry={max_asymmetry:.8f}."
        )

    if not validation_messages:
        validation_messages.append("PASS")

    return row_labels, matrix, validation_messages


def finite_values(values):
    return [v for v in values if not math.isnan(v)]


def finite_mean(values):
    vals = finite_values(values)
    return sum(vals) / len(vals) if vals else math.nan


def finite_min(values):
    vals = finite_values(values)
    return min(vals) if vals else math.nan


def finite_max(values):
    vals = finite_values(values)
    return max(vals) if vals else math.nan


def summarize_segment(segment, path, symmetry_tolerance, diagonal_tolerance):
    labels, identity_matrix, validation_messages = read_identity_matrix(
        path, symmetry_tolerance, diagonal_tolerance
    )
    groups = {g: [i for i, lab in enumerate(labels) if allele(lab) == g] for g in "ABCD"}

    print(f"{segment}: A={len(groups['A'])}, B={len(groups['B'])}, C={len(groups['C'])}, D={len(groups['D'])}")

    validation_rows = [{
        "Segment": segment,
        "Input_file": path.name,
        "N_sequences": len(labels),
        "n_A": len(groups["A"]),
        "n_B": len(groups["B"]),
        "n_C": len(groups["C"]),
        "n_D": len(groups["D"]),
        "Validation_message": msg,
    } for msg in validation_messages]

    B = groups["B"]
    if not B:
        print(f"WARNING: {segment} has no B allele; skipped.")
        return [], validation_rows

    def identity_value(i, j):
        return identity_matrix[i][j]

    def distance_value(i, j):
        v = identity_matrix[i][j]
        return math.nan if math.isnan(v) else 100.0 - v

    def within(indices, func):
        return [func(indices[i], indices[j]) for i in range(len(indices)) for j in range(i + 1, len(indices))]

    def between(indices1, indices2, func):
        return [func(i, j) for i in indices1 for j in indices2]

    B_identity = within(B, identity_value)
    B_distance = within(B, distance_value)
    mean_B_identity = finite_mean(B_identity)
    mean_B_distance = finite_mean(B_distance)

    rows = []
    for other in COMPARISON_ORDER:
        G = groups[other]
        if not G:
            continue

        other_identity = within(G, identity_value)
        other_distance = within(G, distance_value)
        between_identity = between(B, G, identity_value)
        between_distance = between(B, G, distance_value)

        mean_other_identity = finite_mean(other_identity)
        mean_other_distance = finite_mean(other_distance)
        mean_between_identity = finite_mean(between_identity)
        mean_between_distance = finite_mean(between_distance)

        net_between = mean_between_distance - ((mean_B_distance + mean_other_distance) / 2.0)
        identity_ratio = (
            mean_B_identity / mean_between_identity
            if not math.isnan(mean_between_identity) and mean_between_identity != 0
            else math.nan
        )
        distance_ratio = (
            mean_B_distance / mean_between_distance
            if not math.isnan(mean_between_distance) and mean_between_distance != 0
            else math.nan
        )

        rows.append({
            "Segment": segment,
            "Comparison": f"B vs {other}",
            "n_B": len(B),
            "n_other": len(G),
            "Within_B_pair_count": len(finite_values(B_distance)),
            "Within_other_pair_count": len(finite_values(other_distance)),
            "Between_pair_count": len(finite_values(between_distance)),
            "Within_B_identity_%": mean_B_identity,
            "Within_other_identity_%": mean_other_identity,
            "Between_identity_%": mean_between_identity,
            "Within_B_distance_%": mean_B_distance,
            "Within_other_distance_%": mean_other_distance,
            "Between_distance_%": mean_between_distance,
            "Net_between_%": net_between,
            "Identity_ratio": identity_ratio,
            "Distance_ratio": distance_ratio,
            "Within_B_distance_min_%": finite_min(B_distance),
            "Within_B_distance_max_%": finite_max(B_distance),
            "Within_other_distance_min_%": finite_min(other_distance),
            "Within_other_distance_max_%": finite_max(other_distance),
            "Between_distance_min_%": finite_min(between_distance),
            "Between_distance_max_%": finite_max(between_distance),
        })

    return rows, validation_rows


SUMMARY_FIELDS = [
    "Segment", "Comparison", "n_B", "n_other",
    "Within_B_pair_count", "Within_other_pair_count", "Between_pair_count",
    "Within_B_identity_%", "Within_other_identity_%", "Between_identity_%",
    "Within_B_distance_%", "Within_other_distance_%", "Between_distance_%",
    "Net_between_%", "Identity_ratio", "Distance_ratio",
    "Within_B_distance_min_%", "Within_B_distance_max_%",
    "Within_other_distance_min_%", "Within_other_distance_max_%",
    "Between_distance_min_%", "Between_distance_max_%",
]

VALIDATION_FIELDS = [
    "Segment", "Input_file", "N_sequences", "n_A", "n_B", "n_C", "n_D", "Validation_message"
]


def format_value(value):
    if isinstance(value, float):
        return "" if math.isnan(value) else f"{value:.6f}"
    return value


def write_csv(path, rows, fields):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row.get(field, "")) for field in fields})


def make_distance_plot(rows, output_dir):
    if not rows:
        return
    labels = [f"{r['Segment']}\n{r['Comparison']}" for r in rows]
    within_B = np.array([float(r["Within_B_distance_%"]) for r in rows])
    within_other = np.array([float(r["Within_other_distance_%"]) for r in rows])
    between = np.array([float(r["Between_distance_%"]) for r in rows])
    net_between = np.array([float(r["Net_between_%"]) for r in rows])
    x = np.arange(len(labels)); width = 0.20
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.bar(x - 1.5*width, within_B, width, label="Within B")
    ax.bar(x - 0.5*width, within_other, width, label="Within other allele")
    ax.bar(x + 0.5*width, between, width, label="Between groups")
    ax.bar(x + 1.5*width, net_between, width, label="Net between")
    ax.set_ylabel("Genetic distance (%)")
    ax.set_xlabel("Segment and allele comparison")
    ax.set_title("Within- and Between-Allele Genetic Distances")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.legend(); plt.tight_layout()
    plt.savefig(output_dir / "P1_allele_genetic_distance_plot.png", dpi=300, bbox_inches="tight")
    plt.savefig(output_dir / "P1_allele_genetic_distance_plot.pdf", bbox_inches="tight")
    plt.close()


def make_ratio_plot(rows, output_dir):
    if not rows:
        return
    labels = [f"{r['Segment']}\n{r['Comparison']}" for r in rows]
    identity_ratios = np.array([float(r["Identity_ratio"]) for r in rows])
    distance_ratios = np.array([float(r["Distance_ratio"]) for r in rows])
    x = np.arange(len(labels)); width = 0.35
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.bar(x - width/2, identity_ratios, width, label="Identity ratio")
    ax.bar(x + width/2, distance_ratios, width, label="Distance ratio")
    ax.axhline(y=1.0, linestyle="--", linewidth=1)
    ax.set_ylabel("Ratio"); ax.set_xlabel("Segment and allele comparison")
    ax.set_title("Within-B / Between-Allele Identity and Distance Ratios")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.legend(); plt.tight_layout()
    plt.savefig(output_dir / "P1_allele_ratio_plot.png", dpi=300, bbox_inches="tight")
    plt.savefig(output_dir / "P1_allele_ratio_plot.pdf", bbox_inches="tight")
    plt.close()


def make_net_distance_plot(rows, output_dir):
    if not rows:
        return
    labels = [f"{r['Segment']}\n{r['Comparison']}" for r in rows]
    values = np.array([float(r["Net_between_%"]) for r in rows])
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(12, 6))
    bars = ax.bar(x, values)
    ax.set_ylabel("Net between-group genetic distance (%)")
    ax.set_xlabel("Segment and allele comparison")
    ax.set_title("Net Genetic Divergence Between Allele B and Other Alleles")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right")
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(), f"{value:.2f}", ha="center", va="bottom", fontsize=8)
    plt.tight_layout()
    plt.savefig(output_dir / "P1_net_between_distance_plot.png", dpi=300, bbox_inches="tight")
    plt.savefig(output_dir / "P1_net_between_distance_plot.pdf", bbox_inches="tight")
    plt.close()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    all_rows = []
    validation_rows = []
    input_files_used = []

    for segment in SEGMENTS:
        path = args.input_dir / f"{segment}.csv"
        if not path.exists():
            print(f"Not found: {path.name} — skipped.")
            continue
        input_files_used.append(path)
        rows, validations = summarize_segment(
            segment, path, args.symmetry_tolerance, args.diagonal_tolerance
        )
        all_rows.extend(rows)
        validation_rows.extend(validations)

    segment_order = {seg: i for i, seg in enumerate(SEGMENTS)}
    comparison_order = {"B vs A": 0, "B vs C": 1, "B vs D": 2}
    all_rows.sort(key=lambda r: (segment_order[r["Segment"]], comparison_order[r["Comparison"]]))
    validation_rows.sort(key=lambda r: (segment_order[r["Segment"]], r["Validation_message"]))

    write_csv(args.output_dir / "P1_allele_genetic_distance_summary.csv", all_rows, SUMMARY_FIELDS)
    write_csv(args.output_dir / "P1_allele_matrix_validation.csv", validation_rows, VALIDATION_FIELDS)
    write_metadata(args.output_dir, input_files_used, args)

    make_distance_plot(all_rows, args.output_dir)
    make_ratio_plot(all_rows, args.output_dir)
    make_net_distance_plot(all_rows, args.output_dir)

    print("\n" + "="*72)
    print("ANALYSIS COMPLETE")
    print("="*72)
    print(f"Script version:        {SCRIPT_VERSION}")
    print(f"Input files analyzed:  {len(input_files_used)}")
    print(f"Comparisons produced:  {len(all_rows)}")
    print(f"Output directory:      {args.output_dir}")
    print("\nOutput files:")
    print("  P1_allele_genetic_distance_summary.csv")
    print("  P1_allele_matrix_validation.csv")
    print("  P1_distance_run_metadata.txt")
    print("  P1_allele_genetic_distance_plot.png/.pdf")
    print("  P1_allele_ratio_plot.png/.pdf")
    print("  P1_net_between_distance_plot.png/.pdf")


if __name__ == "__main__":
    main()
