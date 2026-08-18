"""Single-person descriptive statistics derived from observed result rows."""

from __future__ import annotations

import json
import math
import re
from collections import Counter
from typing import Any, Iterable

import pandas as pd

BASES = ("A", "C", "G", "T")
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
MISSING_GENOTYPES = {"", ".", "./.", ".|."}
REFERENCE_GENOTYPES = {"0/0", "0|0"}


def _column(frame: pd.DataFrame, *names: str) -> str | None:
    lookup = {str(column).casefold(): str(column) for column in frame.columns}
    for name in names:
        match = lookup.get(name.casefold())
        if match:
            return match
    return None


def _text(value: Any) -> str:
    if value is None or value is pd.NA:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _number(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _round(value: float | None, digits: int = 6) -> float | None:
    return round(value, digits) if value is not None and math.isfinite(value) else None


def _parse_region(value: Any) -> tuple[str, int, int] | None:
    text = _text(value).replace(",", "")
    match = re.fullmatch(r"(?:chr)?([^:\s]+):(\d+)-(\d+)", text, flags=re.IGNORECASE)
    if not match:
        return None
    start, end = int(match.group(2)), int(match.group(3))
    return match.group(1).upper(), min(start, end), max(start, end)


def _region_length(value: Any) -> int | None:
    parsed = _parse_region(value)
    return parsed[2] - parsed[1] + 1 if parsed else None


def _in_region(chromosome: Any, position: Any, region: Any) -> bool:
    parsed = _parse_region(region)
    pos = _number(position)
    if not parsed or pos is None:
        return False
    chrom = _text(chromosome).removeprefix("chr").upper()
    return chrom == parsed[0] and parsed[1] <= int(pos) <= parsed[2]


def _region_category(
    chromosome: Any,
    position: Any,
    *,
    scope_regions: dict[str, Any],
    analysis_scope: str,
    active_region: str,
) -> str:
    pos = _number(position)
    chrom = _text(chromosome)
    if pos is None or not chrom:
        return "unclassified"
    promoter_region = scope_regions.get("promoter_only")
    gene_region = scope_regions.get("gene_only")
    in_promoter = _in_region(chrom, pos, promoter_region)
    in_gene = _in_region(chrom, pos, gene_region)
    if in_promoter and in_gene:
        return "promoter_and_gene"
    if in_promoter:
        return "promoter"
    if in_gene:
        return "gene_body"
    if not promoter_region and analysis_scope == "promoter_only" and _in_region(chrom, pos, active_region):
        return "promoter"
    if not gene_region and analysis_scope == "gene_only" and _in_region(chrom, pos, active_region):
        return "gene_body"
    if _in_region(chrom, pos, active_region):
        return "other_analyzed_region"
    return "unclassified"


def describe_numeric(
    values: Iterable[Any],
    *,
    row_count: int | None = None,
    valid_range: tuple[float, float] | None = None,
) -> dict[str, Any]:
    """Return deterministic descriptive statistics while exposing unusable values."""
    materialized = list(values)
    total = len(materialized) if row_count is None else int(row_count)
    valid: list[float] = []
    missing = 0
    invalid = 0
    for value in materialized:
        text = _text(value)
        if not text:
            missing += 1
            continue
        number = _number(value)
        if number is None or (valid_range and not (valid_range[0] <= number <= valid_range[1])):
            invalid += 1
            continue
        valid.append(number)
    missing += max(0, total - len(materialized))
    result: dict[str, Any] = {
        "row_count": total,
        "valid_count": len(valid),
        "missing_count": missing,
        "invalid_count": invalid,
        "minimum": None,
        "p10": None,
        "q1": None,
        "median": None,
        "mean": None,
        "q3": None,
        "p90": None,
        "maximum": None,
        "iqr": None,
        "population_std": None,
    }
    if not valid:
        result["availability"] = "unavailable"
        result["explanation"] = (
            "The field has no valid numeric values."
            if materialized
            else "The source field is unavailable."
        )
        return result
    series = pd.Series(valid, dtype="float64")
    q1 = float(series.quantile(0.25, interpolation="linear"))
    q3 = float(series.quantile(0.75, interpolation="linear"))
    result.update(
        {
            "minimum": _round(float(series.min())),
            "p10": _round(float(series.quantile(0.10, interpolation="linear"))),
            "q1": _round(q1),
            "median": _round(float(series.median())),
            "mean": _round(float(series.mean())),
            "q3": _round(q3),
            "p90": _round(float(series.quantile(0.90, interpolation="linear"))),
            "maximum": _round(float(series.max())),
            "iqr": _round(q3 - q1),
            "population_std": _round(float(series.std(ddof=0))),
            "availability": "available",
            "explanation": None,
        }
    )
    return result


def _distribution(counter: Counter[str], *, denominator: int) -> list[dict[str, Any]]:
    return [
        {
            "category": category,
            "count": int(count),
            "percentage": _round(100.0 * count / denominator, 3) if denominator else None,
            "denominator": int(denominator),
        }
        for category, count in counter.items()
    ]


def _json_identity(value: Any) -> str:
    try:
        return json.dumps(value, sort_keys=True, default=str, ensure_ascii=False)
    except (TypeError, ValueError):
        return repr(value)


def field_coverage(frame: pd.DataFrame) -> list[dict[str, Any]]:
    """Summarize every raw field without disclosing its distinct values."""
    if frame is None:
        return []
    rows: list[dict[str, Any]] = []
    for column in frame.columns:
        values = frame[column].tolist()
        populated_values = [value for value in values if _text(value)]
        scalar_values = [value for value in populated_values if not isinstance(value, (list, tuple, set, dict))]
        numeric = describe_numeric(scalar_values, row_count=len(values))
        numeric_meaningful = bool(scalar_values) and numeric["valid_count"] == len(scalar_values)
        rows.append(
            {
                "field": str(column),
                "dtype": str(frame[column].dtype),
                "row_count": int(len(values)),
                "populated_count": int(len(populated_values)),
                "missing_count": int(len(values) - len(populated_values)),
                "distinct_count": int(len({_json_identity(value) for value in populated_values})),
                "numeric_summary": numeric if numeric_meaningful else None,
                "numeric_summary_explanation": None
                if numeric_meaningful
                else "The field is non-numeric, mixed-type, or unavailable.",
            }
        )
    return rows


def _genotype_category(row: dict[str, Any]) -> str:
    zygosity = _text(row.get("zygosity")).casefold()
    if zygosity in {"heterozygous", "het"}:
        return "heterozygous"
    if zygosity in {"homozygous_alt", "homozygous alternate", "hom_alt", "homozygous-alternate"}:
        return "homozygous_alternate"
    genotype = _text(row.get("GT") or row.get("gt_raw") or row.get("genotype"))
    if genotype in MISSING_GENOTYPES:
        return "missing"
    if genotype in REFERENCE_GENOTYPES or zygosity in {"homozygous_ref", "reference"}:
        return "reference"
    alleles = re.split(r"[/|]", genotype)
    if len(alleles) == 2 and all(allele.isdigit() for allele in alleles):
        if alleles[0] != alleles[1]:
            return "heterozygous"
        if alleles[0] != "0":
            return "homozygous_alternate"
    return "other"


def _alt_alleles(row: dict[str, Any]) -> list[str]:
    raw = row.get("alt_alleles")
    if isinstance(raw, (list, tuple)):
        return [_text(value).upper() for value in raw if _text(value)]
    text = _text(row.get("ALT") or row.get("alt"))
    return [value.strip().upper() for value in text.split(",") if value.strip() and value.strip() != "."]


def _variant_type(row: dict[str, Any]) -> str:
    ref = _text(row.get("REF") or row.get("ref")).upper()
    alts = _alt_alleles(row)
    if len(alts) > 1:
        return "multiallelic"
    if not ref or not alts:
        return "complex"
    alt = alts[0]
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    if len(ref) < len(alt):
        return "insertion"
    if len(ref) > len(alt):
        return "deletion"
    if len(ref) > 1 and len(alt) > 1:
        return "mnv"
    return "complex"


def _allele_dosages(row: dict[str, Any], alts: list[str]) -> list[float | None]:
    raw = row.get("allele_dosage_per_alt")
    if isinstance(raw, dict):
        values = [raw.get(alt, raw.get(str(index + 1))) for index, alt in enumerate(alts)]
    elif isinstance(raw, (list, tuple)):
        values = list(raw)
    elif _text(raw):
        try:
            parsed = json.loads(_text(raw))
            values = list(parsed.values()) if isinstance(parsed, dict) else list(parsed) if isinstance(parsed, list) else [parsed]
        except (json.JSONDecodeError, TypeError):
            values = [_text(raw)]
    else:
        values = []
    if values:
        return [_number(value) for value in values[: len(alts)]] + [None] * max(0, len(alts) - len(values))
    gt_codes = row.get("gt_codes")
    if isinstance(gt_codes, (list, tuple)):
        return [float(sum(1 for code in gt_codes if _number(code) == index + 1)) for index in range(len(alts))]
    genotype = _text(row.get("GT") or row.get("gt_raw") or row.get("genotype"))
    alleles = re.split(r"[/|]", genotype)
    if alleles and all(allele.isdigit() for allele in alleles):
        return [float(sum(1 for allele in alleles if allele == str(index + 1))) for index in range(len(alts))]
    return [None] * len(alts)


def _field_values(rows: list[dict[str, Any]], names: tuple[str, ...], *, list_first: bool = False) -> list[Any]:
    values: list[Any] = []
    for row in rows:
        lookup = {str(key).casefold(): value for key, value in row.items()}
        value = next((lookup[name.casefold()] for name in names if name.casefold() in lookup), None)
        if list_first and isinstance(value, (list, tuple)):
            value = value[0] if value else None
        values.append(value)
    return values


def build_variant_statistics(
    all_rows: list[dict[str, Any]],
    primary_rows: list[dict[str, Any]],
    *,
    scope_regions: dict[str, Any],
    analysis_scope: str,
    active_region: str,
    raw_frame: pd.DataFrame,
) -> dict[str, Any]:
    gt_categories = Counter(_genotype_category(row) for row in all_rows)
    filter_pass_count = sum(
        1
        for row in all_rows
        if bool(row.get("filter_pass"))
        or _text(row.get("FILTER") or row.get("filter") or row.get("filter_status")).upper() in {"PASS", "."}
    )
    identifiers = [_text(row.get("ID") or row.get("id") or row.get("rsid")) for row in all_rows]
    named_rsid_count = sum(bool(re.fullmatch(r"rs\d+", identifier, flags=re.IGNORECASE)) for identifier in identifiers)
    unlabeled_count = sum(identifier in {"", ".", "Unlabeled in source VCF"} for identifier in identifiers)
    region_counter: Counter[str] = Counter()
    type_counter: Counter[str] = Counter()
    ref_counter: Counter[str] = Counter({base: 0 for base in BASES})
    alt_counter: Counter[str] = Counter({base: 0 for base in BASES})
    dosage_counter: Counter[str] = Counter({base: 0 for base in BASES})
    substitutions: Counter[str] = Counter({f"{ref}>{alt}": 0 for ref in BASES for alt in BASES if ref != alt})
    dosage_evaluable = 0
    transition_count = 0
    transversion_count = 0
    for row in primary_rows:
        chrom = row.get("CHROM") or row.get("chrom") or row.get("chromosome")
        pos = row.get("POS") or row.get("pos") or row.get("position")
        region_counter[
            _region_category(
                chrom,
                pos,
                scope_regions=scope_regions,
                analysis_scope=analysis_scope,
                active_region=active_region,
            )
        ] += 1
        variant_type = _variant_type(row)
        type_counter[variant_type] += 1
        ref = _text(row.get("REF") or row.get("ref")).upper()
        alts = _alt_alleles(row)
        if len(ref) != 1 or ref not in BASES or not alts or any(len(alt) != 1 or alt not in BASES for alt in alts):
            continue
        ref_counter[ref] += 1
        dosages = _allele_dosages(row, alts)
        if any(dosage is not None for dosage in dosages):
            dosage_evaluable += 1
        for alt, dosage in zip(alts, dosages):
            alt_counter[alt] += 1
            substitutions[f"{ref}>{alt}"] += 1
            if (ref, alt) in TRANSITIONS:
                transition_count += 1
            else:
                transversion_count += 1
            if dosage is not None:
                dosage_counter[alt] += dosage

    region_order = ("promoter", "gene_body", "promoter_and_gene", "other_analyzed_region", "unclassified")
    type_order = ("snv", "insertion", "deletion", "mnv", "complex", "multiallelic")
    genotype_order = ("heterozygous", "homozygous_alternate", "other", "reference", "missing")
    region_rows = _distribution(Counter({key: region_counter[key] for key in region_order}), denominator=len(primary_rows))
    density: list[dict[str, Any]] = []
    for key, label, count in (
        ("promoter_only", "promoter", region_counter["promoter"] + region_counter["promoter_and_gene"]),
        ("gene_only", "gene_body", region_counter["gene_body"] + region_counter["promoter_and_gene"]),
        ("active", "active_region", len(primary_rows)),
    ):
        region = active_region if key == "active" else scope_regions.get(key)
        length = _region_length(region)
        density.append(
            {
                "region": label,
                "coordinates": _text(region),
                "length_bp": length,
                "variant_count": int(count),
                "variants_per_kb": _round(count / (length / 1000.0), 6) if length else None,
                "explanation": None if length else "Region coordinates are unavailable.",
            }
        )

    quality_specs = (
        ("qual", ("QUAL", "qual"), False),
        ("genotype_quality", ("GQ", "gq", "genotype_quality"), False),
        ("depth", ("DP", "dp", "depth"), False),
        ("confidence_score", ("confidence_score",), False),
        ("alternate_allele_fraction", ("sample_af", "AF", "allele_fraction"), True),
    )
    quality = [
        {"metric": label, **describe_numeric(_field_values(primary_rows, names, list_first=list_first), row_count=len(primary_rows))}
        for label, names, list_first in quality_specs
    ]
    counts = {
        "total_count": int(len(all_rows)),
        "raw_row_count": int(len(all_rows)),
        "filter_pass_count": int(filter_pass_count),
        "filter_non_pass_count": int(len(all_rows) - filter_pass_count),
        "qc_passing_non_reference_count": int(len(primary_rows)),
        "non_reference_count": int(sum(bool(row.get("non_reference")) for row in all_rows)),
        "reference_count": int(gt_categories["reference"]),
        "missing_genotype_count": int(gt_categories["missing"]),
        "named_rsid_count": int(named_rsid_count),
        "named_variant_count": int(named_rsid_count),
        "unlabeled_count": int(unlabeled_count),
        "unlabeled_variant_count": int(unlabeled_count),
    }
    return {
        "status": "descriptive" if all_rows else "no_data",
        "counts": counts,
        "by_region": region_rows,
        "by_type": _distribution(Counter({key: type_counter[key] for key in type_order}), denominator=len(primary_rows)),
        "genotypes": _distribution(Counter({key: gt_categories[key] for key in genotype_order}), denominator=len(all_rows)),
        "reference_bases": _distribution(ref_counter, denominator=sum(ref_counter.values())),
        "alternate_bases": _distribution(alt_counter, denominator=sum(alt_counter.values())),
        "dosage_weighted_alternate_copies": [
            {"base": base, "allele_copies": _round(float(dosage_counter[base]), 6)} for base in BASES
        ],
        "dosage_evaluable_snv_count": int(dosage_evaluable),
        "substitutions": _distribution(substitutions, denominator=sum(substitutions.values())),
        "transition_transversion": {
            "transition_count": int(transition_count),
            "transversion_count": int(transversion_count),
            "ratio": _round(transition_count / transversion_count, 6) if transversion_count else None,
        },
        "density": density,
        "quality": quality,
        "raw_field_coverage": field_coverage(raw_frame),
    }


def _gene_named_mask(frame: pd.DataFrame, gene: str) -> pd.Series:
    mask = pd.Series(False, index=frame.index, dtype="bool")
    normalized = gene.strip().upper()
    if not normalized:
        return mask
    for column_name in ("GencodeBasicV12_NAME", "UCSC_RefGene_Name", "gene"):
        column = _column(frame, column_name)
        if not column:
            continue
        column_mask = frame[column].apply(
            lambda value: normalized in {token.strip().upper() for token in re.split(r"[;,]", _text(value)) if token.strip()}
        )
        mask |= column_mask
    return mask


def _beta_values(frame: pd.DataFrame, beta_column: str | None) -> list[Any]:
    return frame[beta_column].tolist() if beta_column and beta_column in frame else [None] * len(frame)


def _subset_summary(label: str, frame: pd.DataFrame, beta_column: str | None) -> dict[str, Any]:
    return {"subset": label, **describe_numeric(_beta_values(frame, beta_column), row_count=len(frame), valid_range=(0.0, 1.0))}


def _categorical_beta_groups(frame: pd.DataFrame, column: str | None, beta_column: str | None) -> list[dict[str, Any]]:
    if not column:
        return []
    groups: dict[str, list[int]] = {}
    for index, value in frame[column].items():
        tokens = [token.strip() for token in re.split(r"[;]", _text(value)) if token.strip()] or ["Unannotated"]
        for token in dict.fromkeys(tokens):
            groups.setdefault(token, []).append(index)
    return [
        {"category": category, **_subset_summary(category, frame.loc[indices], beta_column)}
        for category, indices in sorted(groups.items(), key=lambda item: (-len(item[1]), item[0]))
    ]


def _beta_histogram(values: Iterable[Any]) -> list[dict[str, Any]]:
    valid = [number for value in values if (number := _number(value)) is not None and 0.0 <= number <= 1.0]
    counts = [0] * 10
    for value in valid:
        counts[min(9, int(value * 10))] += 1
    return [
        {
            "bin": f"{index / 10:.1f}-{(index + 1) / 10:.1f}",
            "lower_inclusive": index / 10,
            "upper_inclusive": index == 9,
            "count": int(count),
            "percentage": _round(100.0 * count / len(valid), 3) if valid else None,
            "denominator": len(valid),
        }
        for index, count in enumerate(counts)
    ]


def _annotation_coverage(frame: pd.DataFrame) -> list[dict[str, Any]]:
    specifications = (
        ("gene_name", ("GencodeBasicV12_NAME", "UCSC_RefGene_Name", "gene")),
        ("refgene_group", ("UCSC_RefGene_Group",)),
        ("cpg_island", ("UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island")),
        ("enhancer", ("Phantom4_Enhancers", "enhancer")),
        ("dnase_hypersensitivity", ("DNase_Hypersensitivity_NAME", "DNase_Hypersensitivity_Evidence_Count")),
        ("regulatory_feature", ("Regulatory_Feature_Name", "Regulatory_Feature_Group", "OpenChromatin_NAME")),
        ("transcription_factor_binding", ("TFBS_NAME", "TFBS_Evidence_Count")),
    )
    rows: list[dict[str, Any]] = []
    for label, names in specifications:
        columns = [column for name in names if (column := _column(frame, name))]
        if not columns:
            rows.append(
                {
                    "annotation": label,
                    "available": False,
                    "annotated_count": None,
                    "percentage": None,
                    "explanation": "No supported source column is available.",
                }
            )
            continue
        mask = pd.Series(False, index=frame.index, dtype="bool")
        for column in columns:
            mask |= frame[column].apply(lambda value: bool(_text(value)) and _text(value).casefold() not in {"none", "nan", "."})
        count = int(mask.sum())
        rows.append(
            {
                "annotation": label,
                "available": True,
                "columns": columns,
                "annotated_count": count,
                "percentage": _round(100.0 * count / len(frame), 3) if len(frame) else None,
                "explanation": None,
            }
        )
    return rows


def _probe_extremes(frame: pd.DataFrame, beta_column: str | None, *, limit: int = 5) -> dict[str, Any]:
    if not beta_column or frame.empty:
        return {"subset": "qc_included", "lowest": [], "highest": []}
    rows: list[dict[str, Any]] = []
    probe_column = _column(frame, "probe_id", "IlmnID", "Name")
    chrom_column = _column(frame, "chrom", "CHR")
    pos_column = _column(frame, "pos", "MAPINFO")
    gene_column = _column(frame, "GencodeBasicV12_NAME", "UCSC_RefGene_Name", "gene")
    group_column = _column(frame, "UCSC_RefGene_Group")
    for _, row in frame.iterrows():
        beta = _number(row.get(beta_column))
        if beta is None or not 0.0 <= beta <= 1.0:
            continue
        rows.append(
            {
                "probe_id": _text(row.get(probe_column)) if probe_column else "",
                "beta": _round(beta),
                "chromosome": _text(row.get(chrom_column)) if chrom_column else "",
                "position": int(_number(row.get(pos_column))) if pos_column and _number(row.get(pos_column)) is not None else None,
                "gene_annotation": _text(row.get(gene_column)) if gene_column else "",
                "gene_region_annotation": _text(row.get(group_column)) if group_column else "",
                "region": _text(row.get("single_person_region")),
            }
        )
    ordered = sorted(rows, key=lambda item: (item["beta"], item["probe_id"]))
    return {"subset": "qc_included", "lowest": ordered[:limit], "highest": list(reversed(ordered[-limit:]))}


def build_methylation_statistics(
    all_frame: pd.DataFrame,
    primary_frame: pd.DataFrame,
    *,
    raw_frame: pd.DataFrame,
    gene: str,
    curated_probe_ids: Iterable[str],
    scope_regions: dict[str, Any],
    analysis_scope: str,
    active_region: str,
) -> dict[str, Any]:
    all_frame = all_frame.copy()
    primary_frame = primary_frame.copy()
    beta_column = _column(all_frame, "beta_value", "beta", "Beta_value")
    primary_beta_column = _column(primary_frame, "beta_value", "beta", "Beta_value")
    chrom_column = _column(all_frame, "chrom", "CHR")
    pos_column = _column(all_frame, "pos", "MAPINFO")
    categories: list[str] = []
    for _, row in all_frame.iterrows():
        categories.append(
            _region_category(
                row.get(chrom_column) if chrom_column else None,
                row.get(pos_column) if pos_column else None,
                scope_regions=scope_regions,
                analysis_scope=analysis_scope,
                active_region=active_region,
            )
        )
    all_frame["single_person_region"] = categories
    if not primary_frame.empty:
        primary_chrom = _column(primary_frame, "chrom", "CHR")
        primary_pos = _column(primary_frame, "pos", "MAPINFO")
        primary_frame["single_person_region"] = [
            _region_category(
                row.get(primary_chrom) if primary_chrom else None,
                row.get(primary_pos) if primary_pos else None,
                scope_regions=scope_regions,
                analysis_scope=analysis_scope,
                active_region=active_region,
            )
            for _, row in primary_frame.iterrows()
        ]

    gene_named = all_frame[_gene_named_mask(all_frame, gene)]
    probe_column = _column(all_frame, "probe_id", "IlmnID", "Name")
    whitelist = {str(probe).strip() for probe in curated_probe_ids if str(probe).strip()}
    curated = all_frame[all_frame[probe_column].astype(str).isin(whitelist)] if probe_column and whitelist else all_frame.iloc[0:0]
    promoter = all_frame[all_frame["single_person_region"].isin(["promoter", "promoter_and_gene"])]
    gene_body = all_frame[all_frame["single_person_region"].isin(["gene_body", "promoter_and_gene"])]
    subsets = [
        _subset_summary("all_rows", all_frame, beta_column),
        _subset_summary("qc_included", primary_frame, primary_beta_column),
        _subset_summary("gene_named", gene_named, beta_column),
        _subset_summary("curated_whitelist", curated, beta_column),
        _subset_summary("promoter", promoter, beta_column),
        _subset_summary("gene_body", gene_body, beta_column),
    ]
    region_groups = [
        {"category": category, **_subset_summary(category, group, beta_column)}
        for category, group in all_frame.groupby("single_person_region", dropna=False, sort=True)
    ]
    quality: list[dict[str, Any]] = []
    for label, names, valid_range in (
        ("m_value", ("m_value", "M_value"), None),
        ("detection_p", ("detection_p", "Detection Pval", "detection_p_value"), (0.0, 1.0)),
        ("bead_count", ("bead_count", "NBeads", "beads"), (0.0, float("inf"))),
    ):
        column = _column(all_frame, *names)
        quality.append(
            {
                "metric": label,
                "source_column": column,
                **describe_numeric(all_frame[column].tolist() if column else [], row_count=len(all_frame), valid_range=valid_range),
            }
        )
    all_beta = describe_numeric(_beta_values(all_frame, beta_column), row_count=len(all_frame), valid_range=(0.0, 1.0))
    return {
        "status": "descriptive" if len(all_frame) else "no_data",
        "counts": {
            "raw_probe_count": int(len(all_frame)),
            "qc_included_probe_count": int(len(primary_frame)),
            "gene_named_probe_count": int(len(gene_named)),
            "curated_whitelist_probe_count": int(len(curated)),
            "promoter_probe_count": int(len(promoter)),
            "gene_body_probe_count": int(len(gene_body)),
            "valid_beta_count": int(all_beta["valid_count"]),
            "missing_beta_count": int(all_beta["missing_count"]),
            "invalid_beta_count": int(all_beta["invalid_count"]),
        },
        "subsets": subsets,
        "by_region": region_groups,
        "by_refgene_group": _categorical_beta_groups(
            all_frame, _column(all_frame, "UCSC_RefGene_Group"), beta_column
        ),
        "by_cpg_island_relation": _categorical_beta_groups(
            all_frame, _column(all_frame, "Relation_to_UCSC_CpG_Island"), beta_column
        ),
        "beta_histogram": _beta_histogram(_beta_values(all_frame, beta_column)),
        "quality": quality,
        "annotation_coverage": _annotation_coverage(all_frame),
        "extremes": _probe_extremes(all_frame, beta_column),
        "raw_field_coverage": field_coverage(raw_frame),
    }


def build_single_person_statistics(
    *,
    gene: str,
    variants: pd.DataFrame,
    all_variant_rows: list[dict[str, Any]],
    primary_variant_rows: list[dict[str, Any]],
    methylation: pd.DataFrame,
    all_methylation_rows: list[dict[str, Any]],
    primary_methylation_rows: list[dict[str, Any]],
    scope_regions: dict[str, Any] | None,
    analysis_scope: str,
    active_region: str,
    curated_probe_ids: Iterable[str] = (),
) -> dict[str, Any]:
    """Build the standard report's statistics section from one person's rows."""
    regions = dict(scope_regions or {})
    all_methylation = pd.DataFrame(all_methylation_rows) if all_methylation_rows else methylation.iloc[0:0].copy()
    primary_methylation = (
        pd.DataFrame(primary_methylation_rows) if primary_methylation_rows else methylation.iloc[0:0].copy()
    )
    variant_statistics = build_variant_statistics(
        all_variant_rows,
        primary_variant_rows,
        scope_regions=regions,
        analysis_scope=analysis_scope,
        active_region=active_region,
        raw_frame=variants,
    )
    methylation_statistics = build_methylation_statistics(
        all_methylation,
        primary_methylation,
        raw_frame=methylation,
        gene=gene,
        curated_probe_ids=curated_probe_ids,
        scope_regions=regions,
        analysis_scope=analysis_scope,
        active_region=active_region,
    )
    has_variants = bool(all_variant_rows)
    has_methylation = bool(all_methylation_rows)
    status = "descriptive" if has_variants and has_methylation else "partial" if has_variants or has_methylation else "no_data"
    summary_record = {
        "entity_type": "gene",
        "entity_key": gene,
        "family": "single_sample_descriptive",
        "method": "raw_observation_summary",
        "status": status,
        "variant_statistics": variant_statistics,
        "methylation_statistics": methylation_statistics,
        "limitations": [
            "Descriptive statistics for one person; no population comparison or significance test.",
            "Unavailable input fields remain null and are not imputed.",
        ],
    }
    return {
        "scope": "single_person",
        "status": status,
        "gene": gene,
        "variant_statistics": variant_statistics,
        "methylation_statistics": methylation_statistics,
        "records": [summary_record],
        "limitations": summary_record["limitations"],
    }
