"""Matched-reference exploratory statistics with explicit safeguards."""

from __future__ import annotations

import math
from collections import defaultdict
from statistics import median
from typing import Any, Iterable

REFERENCE_MINIMUM = 30
MATCH_FIELDS = ("tissue", "platform", "normalization", "genome_build")


def reference_compatibility(sample_context: dict[str, Any], cohort_context: dict[str, Any]) -> dict[str, Any]:
    mismatches: list[str] = []
    missing: list[str] = []
    for field in MATCH_FIELDS:
        sample_value = str(sample_context.get(field) or "").strip().casefold()
        cohort_value = str(cohort_context.get(field) or "").strip().casefold()
        if not sample_value or not cohort_value:
            missing.append(field)
        elif sample_value != cohort_value:
            mismatches.append(field)
    return {
        "compatible": not mismatches and not missing,
        "mismatches": mismatches,
        "missing_context": missing,
        "policy": "exact tissue/platform/normalization/build match; cohorts are never silently pooled",
    }


def benjamini_hochberg(p_values: Iterable[float | None]) -> list[float | None]:
    values = list(p_values)
    indexed = sorted(
        ((index, float(value)) for index, value in enumerate(values) if value is not None and math.isfinite(float(value))),
        key=lambda item: item[1],
    )
    adjusted: list[float | None] = [None] * len(values)
    running = 1.0
    total = len(indexed)
    for reverse_rank, (index, value) in enumerate(reversed(indexed), start=1):
        rank = total - reverse_rank + 1
        running = min(running, value * total / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def empirical_single_sample_result(value: float, reference_values: Iterable[float]) -> dict[str, Any]:
    reference = sorted(float(item) for item in reference_values if item is not None and math.isfinite(float(item)))
    if len(reference) < REFERENCE_MINIMUM:
        return {
            "status": "descriptive_only",
            "reference_sample_count": len(reference),
            "minimum_required": REFERENCE_MINIMUM,
            "limitations": ["reference_cohort_below_30"],
        }
    center = median(reference)
    deviations = [abs(item - center) for item in reference]
    mad = median(deviations)
    robust_effect = (value - center) / (1.4826 * mad) if mad > 0 else None
    below = sum(item <= value for item in reference)
    above = sum(item >= value for item in reference)
    percentile = 100.0 * below / len(reference)
    raw_p = min(1.0, 2.0 * min((below + 1) / (len(reference) + 1), (above + 1) / (len(reference) + 1)))
    return {
        "status": "exploratory",
        "method": "empirical_two_sided_tail",
        "value": value,
        "reference_sample_count": len(reference),
        "reference_median": center,
        "reference_mad": mad,
        "effect_size": robust_effect,
        "percentile": percentile,
        "raw_p": raw_p,
        "limitations": ["hypothesis_generating_single_sample_comparison"],
    }


def apply_family_fdr(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        grouped[str(row.get("family") or "unspecified")].append(index)
    result = [dict(row) for row in rows]
    for family_indices in grouped.values():
        q_values = benjamini_hochberg(result[index].get("raw_p") for index in family_indices)
        for index, q_value in zip(family_indices, q_values):
            result[index]["q_value"] = q_value
    return result
