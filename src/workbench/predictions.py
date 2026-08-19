"""Single-person prediction preparation and source-native annotations.

This module deliberately separates measured/person-matched source annotations
from executable model output.  It contains no network calls.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any, Iterable


SUPPORTED_ALPHAGENOME_MODALITIES = (
    "RNA_SEQ",
    "CAGE",
    "DNASE",
    "ATAC",
    "CHIP_HISTONE",
    "CHIP_TF",
    "SPLICE_SITES",
    "SPLICE_JUNCTIONS",
)
REGULATORY_OVERVIEW_MODALITIES = SUPPORTED_ALPHAGENOME_MODALITIES[:6]
# Intersection of the canonical request helper's accepted lengths and the
# reviewed AlphaGenome 0.8.0 SDK's enumerated lengths.
SUPPORTED_SEQUENCE_LENGTHS = (131_072, 524_288, 1_048_576)
DEFAULT_SEQUENCE_LENGTH = 1_048_576
MAX_MODEL_VARIANTS = 20
MAX_ONTOLOGY_TERMS = 5


def _text(value: Any) -> str:
    return "" if value is None else str(value).strip()


def _number(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def normalize_build(value: Any) -> str:
    normalized = _text(value).casefold().replace("_", "").replace("-", "")
    if normalized in {"hg38", "grch38", "b38"}:
        return "GRCh38"
    if normalized in {"hg19", "grch37", "b37"}:
        return "GRCh37"
    return _text(value)


def normalize_chromosome(value: Any) -> str:
    chromosome = re.sub(r"^chr", "", _text(value), flags=re.IGNORECASE)
    return chromosome.upper() if chromosome.upper() in {"X", "Y", "M", "MT"} else chromosome


def _first(record: dict[str, Any], *names: str) -> Any:
    lookup = {str(key).casefold(): value for key, value in record.items()}
    for name in names:
        if name.casefold() in lookup:
            return lookup[name.casefold()]
    return None


def _interval(value: Any) -> tuple[str, int, int] | None:
    if isinstance(value, dict):
        chromosome = normalize_chromosome(value.get("chromosome") or value.get("chrom"))
        start = _number(value.get("start"))
        end = _number(value.get("end") or value.get("stop"))
        if chromosome and start is not None and end is not None:
            return chromosome, int(start), int(end)
    match = re.fullmatch(r"(?:chr)?([^:\s]+):([\d,]+)-([\d,]+)", _text(value))
    if not match:
        return None
    chromosome, start, end = match.groups()
    return normalize_chromosome(chromosome), int(start.replace(",", "")), int(end.replace(",", ""))


def _region_for_locus(chromosome: str, position: int, scope_regions: dict[str, Any]) -> str:
    promoter = _interval(scope_regions.get("promoter_only") or scope_regions.get("promoter"))
    gene = _interval(scope_regions.get("gene_only") or scope_regions.get("gene_body") or scope_regions.get("gene"))
    in_promoter = bool(promoter and promoter[0] == chromosome and promoter[1] <= position <= promoter[2])
    in_gene = bool(gene and gene[0] == chromosome and gene[1] <= position <= gene[2])
    if in_promoter and in_gene:
        return "promoter_and_gene"
    if in_promoter:
        return "promoter"
    if in_gene:
        return "gene_body"
    return "other_analyzed_region"


def _called_alt_alleles(row: dict[str, Any], alternates: list[str]) -> list[str]:
    dosage = _first(row, "allele_dosage_per_alt")
    if isinstance(dosage, dict):
        called = [str(alt).upper() for alt, count in dosage.items() if (_number(count) or 0) > 0]
        if called:
            return called
    genotype_alleles = _first(row, "genotype_alleles")
    if isinstance(genotype_alleles, list):
        called = {
            _text(allele).upper()
            for allele in genotype_alleles
            if _text(allele) and _text(allele).upper() in alternates
        }
        if called:
            return [alt for alt in alternates if alt in called]
    codes = _first(row, "gt_codes")
    if isinstance(codes, list):
        called_indexes = {int(code) for code in codes if isinstance(code, (int, float)) and int(code) > 0}
        called = [alt for index, alt in enumerate(alternates, start=1) if index in called_indexes]
        if called:
            return called
    raw_genotype = _text(_first(row, "GT", "gt_raw", "genotype"))
    if raw_genotype and raw_genotype not in {".", "./.", ".|."}:
        tokens = re.split(r"[/|]", raw_genotype.split(":", 1)[0])
        called_indexes: set[int] = set()
        for token in tokens:
            if token in {"", "."}:
                continue
            try:
                index = int(token)
            except ValueError:
                return []
            if index > 0:
                called_indexes.add(index)
        called = [alt for index, alt in enumerate(alternates, start=1) if index in called_indexes]
        return called
    return alternates


def observed_variant_alleles(
    rows: Iterable[dict[str, Any]],
    *,
    genome_build: str,
    scope_regions: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Return normalized observed ALT alleles from QC-passing non-reference rows."""
    build = normalize_build(genome_build)
    normalized: dict[str, dict[str, Any]] = {}
    for row in rows:
        if not isinstance(row, dict) or row.get("qc_pass") is not True or row.get("non_reference") is not True:
            continue
        chromosome = normalize_chromosome(_first(row, "CHROM", "chromosome", "chrom"))
        position_value = _number(_first(row, "POS", "position", "pos"))
        reference = _text(_first(row, "REF", "reference", "ref")).upper()
        raw_alt = _first(row, "ALT", "alternate", "alt")
        alternates = [item.strip().upper() for item in _text(raw_alt).split(",") if item.strip()]
        if not chromosome or position_value is None or not reference or not alternates:
            continue
        position = int(position_value)
        rsid = _text(_first(row, "ID", "rsid", "variant_id"))
        if not rsid.casefold().startswith("rs"):
            rsid = ""
        region = _text(_first(row, "region_classification", "region_category")) or _region_for_locus(
            chromosome, position, dict(scope_regions or {})
        )
        curated = bool(rsid) or "knowledge base" in _text(_first(row, "id_source")).casefold()
        for alternate in _called_alt_alleles(row, alternates):
            if alternate == reference:
                continue
            key = f"{build}|{chromosome}|{position}|{reference}|{alternate}"
            normalized[key] = {
                "entity_key": key,
                "assembly": build,
                "chromosome": f"chr{chromosome}",
                "position": position,
                "reference": reference,
                "alternate": alternate,
                "variant": f"chr{chromosome}:{position}:{reference}>{alternate}",
                "rsid": rsid,
                "region": region,
                "curated_or_named": curated,
                "genotype": _text(_first(row, "GT", "gt_raw", "genotype")),
                "gq": _number(_first(row, "GQ", "genotype_quality")),
                "dp": _number(_first(row, "DP", "depth")),
                "confidence": _number(_first(row, "confidence_score", "confidence")),
            }
    return sorted(normalized.values(), key=lambda item: item["entity_key"])


def _record_containers(payload: dict[str, Any]) -> Iterable[tuple[dict[str, Any], str]]:
    seen: set[int] = set()
    for container_name, container in (
        ("run", payload),
        ("curated", payload.get("knowledge_base")),
        ("dynamic", payload.get("dynamic_knowledge_base")),
    ):
        if not isinstance(container, dict):
            continue
        for key in ("source_records", "records", "variant_records", "population_records"):
            values = container.get(key)
            if not isinstance(values, list):
                continue
            for record in values:
                if isinstance(record, dict) and id(record) not in seen:
                    seen.add(id(record))
                    yield record, f"{container_name}_{key}"


def _record_locus(record: dict[str, Any]) -> tuple[str, str, int, str, str] | None:
    build = normalize_build(_first(record, "reference_genome", "genome_build", "assembly", "build"))
    chromosome = normalize_chromosome(_first(record, "chromosome", "chrom"))
    position_value = _number(_first(record, "position", "pos"))
    reference = _text(_first(record, "reference", "ref")).upper()
    alternate = _text(_first(record, "alternate", "alt")).upper()
    if build and chromosome and position_value is not None and reference and alternate:
        return build, chromosome, int(position_value), reference, alternate
    source_id = _text(_first(record, "variant_id", "source_id"))
    match = re.fullmatch(r"(?:chr)?([^:\-]+)[-:](\d+)[-:]([ACGTN]+)[-:]([ACGTN]+)", source_id, re.IGNORECASE)
    if match and build:
        chromosome, position, reference, alternate = match.groups()
        return build, normalize_chromosome(chromosome), int(position), reference.upper(), alternate.upper()
    return None


def _predictor_rows(record: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    values = record.get("in_silico_predictors") or record.get("predictors")
    if isinstance(values, dict):
        values = [{"id": key, "value": value} for key, value in values.items()]
    if isinstance(values, list):
        rows.extend(item for item in values if isinstance(item, dict))
    transcript = record.get("transcript_consequence")
    if isinstance(transcript, dict):
        for key in ("polyphen_prediction", "sift_prediction"):
            if transcript.get(key) not in (None, ""):
                rows.append({"id": key.replace("_prediction", ""), "value": transcript[key]})
    return rows


def mark_exact_source_matches(
    payload: dict[str, Any], observed_alleles: Iterable[dict[str, Any]]
) -> list[dict[str, Any]]:
    """Mark exact source-record matches for deterministic model prioritization."""
    source_matches: dict[str, str] = {}
    for record, _origin in _record_containers(payload):
        locus = _record_locus(record)
        if locus is None:
            continue
        build, chromosome, position, reference, alternate = locus
        key = f"{build}|{chromosome}|{position}|{reference}|{alternate}"
        source_rsid = _text(record.get("rsid") or record.get("label"))
        source_matches[key] = source_rsid if source_rsid.casefold().startswith("rs") else ""
    enriched = []
    for item in observed_alleles:
        row = dict(item)
        if row.get("entity_key") in source_matches:
            row["curated_or_named"] = True
            row["exact_source_match"] = True
            if not row.get("rsid") and source_matches[row["entity_key"]]:
                row["rsid"] = source_matches[row["entity_key"]]
        else:
            row["exact_source_match"] = False
        enriched.append(row)
    return enriched


def source_native_annotations(
    payload: dict[str, Any], observed_alleles: Iterable[dict[str, Any]]
) -> tuple[list[dict[str, Any]], int]:
    """Match source predictor records to observed alleles by exact build and locus."""
    observed = {item["entity_key"]: item for item in observed_alleles}
    annotations: dict[str, dict[str, Any]] = {}
    unmatched_records = 0
    for record, origin in _record_containers(payload):
        predictors = _predictor_rows(record)
        if not predictors:
            continue
        locus = _record_locus(record)
        if not locus:
            unmatched_records += 1
            continue
        build, chromosome, position, reference, alternate = locus
        key = f"{build}|{chromosome}|{position}|{reference}|{alternate}"
        allele = observed.get(key)
        if allele is None:
            unmatched_records += 1
            continue
        for predictor in predictors:
            name = _text(predictor.get("id") or predictor.get("name") or predictor.get("predictor"))
            original = predictor.get("value") if "value" in predictor else predictor.get("score")
            if not name or original in (None, ""):
                continue
            identity_payload = {
                "entity_key": key,
                "predictor": name.casefold(),
                "source": _text(record.get("source_key") or record.get("source")),
                "value": original,
            }
            annotation_id = hashlib.sha256(
                json.dumps(identity_payload, sort_keys=True, default=str).encode("utf-8")
            ).hexdigest()[:24]
            annotations[annotation_id] = {
                "annotation_id": annotation_id,
                "entity_type": "observed_variant",
                "entity_key": key,
                "variant": allele["variant"],
                "rsid": _text(record.get("rsid")) or allele.get("rsid", ""),
                "predictor": name,
                "native_value": _number(original),
                "native_value_original": original,
                "native_label": _text(predictor.get("label")),
                "flags": list(
                    dict.fromkeys(
                        str(flag)
                        for values in (
                            predictor.get("flags") or [],
                            record.get("flags") or [],
                            record.get("filters") or [],
                        )
                        for flag in (values if isinstance(values, (list, tuple, set)) else [values])
                        if str(flag).strip()
                    )
                ),
                "source_key": _text(record.get("source_key") or record.get("source")),
                "source": _text(record.get("source") or record.get("source_key")),
                "source_release": _text(
                    record.get("source_release") or record.get("release") or record.get("dataset")
                ),
                "source_record_id": _text(record.get("record_id") or record.get("source_id")),
                "url": _text(record.get("url")),
                "evidence_origin": origin,
                "original_representation": json.loads(json.dumps(predictor, default=str)),
                "match_basis": "exact_assembly_chromosome_position_ref_alt",
                "limitations": [
                    "Source-native computational annotation; scales and directions are predictor-specific.",
                    "Not a diagnosis, pathogenicity classification, phenotype prediction, or statistical significance measurement.",
                ],
            }
    return sorted(
        annotations.values(),
        key=lambda item: (item["variant"], item["predictor"].casefold(), item["source_key"]),
    ), unmatched_records


def select_alphagenome_variants(
    alleles: Iterable[dict[str, Any]], *, maximum: int = MAX_MODEL_VARIANTS
) -> dict[str, Any]:
    rows = [dict(item) for item in alleles]
    if rows and any(item.get("assembly") != "GRCh38" for item in rows):
        return {
            "status": "blocked",
            "selected": [],
            "omitted": [{**item, "omission_reason": "alphagenome_requires_grch38"} for item in rows],
            "blockers": ["unsupported_build"],
            "maximum": maximum,
            "selection_policy": "qc_nonreference_grch38_curated_region_quality_locus",
        }
    eligible_rows = [
        item
        for item in rows
        if re.fullmatch(r"[ACGTN]+", _text(item.get("reference")).upper())
        and re.fullmatch(r"[ACGTN]+", _text(item.get("alternate")).upper())
        and _text(item.get("reference")).upper() != _text(item.get("alternate")).upper()
    ]
    unsupported = [
        {**item, "omission_reason": "unsupported_variant_representation"}
        for item in rows
        if item not in eligible_rows
    ]
    region_priority = {"promoter_and_gene": 0, "promoter": 1, "gene_body": 2, "other_analyzed_region": 3}
    eligible_rows.sort(
        key=lambda item: (
            0 if item.get("curated_or_named") else 1,
            region_priority.get(_text(item.get("region")), 4),
            -(_number(item.get("gq")) or -1),
            -(_number(item.get("dp")) or -1),
            item.get("entity_key", ""),
        )
    )
    selected = eligible_rows[:maximum]
    omitted = [
        *unsupported,
        *[{**item, "omission_reason": "maximum_variant_limit"} for item in eligible_rows[maximum:]],
    ]
    return {
        "status": "available" if selected else "no_data",
        "selected": selected,
        "omitted": omitted,
        "blockers": [] if selected else [
            "no_alphagenome_compatible_variants" if rows else "no_qc_passing_nonreference_variants"
        ],
        "maximum": maximum,
        "selection_policy": "qc_nonreference_grch38_curated_region_quality_locus",
    }


def validate_reference_alleles(selection: dict[str, Any], fasta_path: Path) -> dict[str, Any]:
    """Independently verify selected REF bases through a standard FASTA index."""
    result = {
        **selection,
        "selected": [dict(item) for item in selection.get("selected", [])],
        "omitted": [dict(item) for item in selection.get("omitted", [])],
    }
    fasta = Path(fasta_path)
    index_path = Path(f"{fasta}.fai")
    if not fasta.is_file() or not index_path.is_file():
        result["status"] = "blocked"
        result["blockers"] = sorted(set([*result.get("blockers", []), "reference_fasta_unavailable"]))
        return result
    index: dict[str, tuple[int, int, int, int]] = {}
    for line in index_path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        name, length, offset, line_bases, line_width = fields[:5]
        index[normalize_chromosome(name)] = (int(length), int(offset), int(line_bases), int(line_width))
    blockers = list(result.get("blockers", []))
    validated_selected: list[dict[str, Any]] = []
    with fasta.open("rb") as handle:
        for item in result["selected"]:
            chromosome = normalize_chromosome(item.get("chromosome"))
            metadata = index.get(chromosome)
            if metadata is None:
                item.update(reference_allele_verified=False, reference_validation_status="chromosome_missing")
                blockers.append(f"reference_chromosome_missing:{chromosome}")
                result["omitted"].append({**item, "omission_reason": "reference_chromosome_missing"})
                continue
            length, offset, line_bases, line_width = metadata
            position = int(item["position"])
            reference = str(item["reference"]).upper()
            observed = []
            for delta in range(len(reference)):
                zero_based = position - 1 + delta
                if zero_based < 0 or zero_based >= length:
                    break
                byte_position = offset + (zero_based // line_bases) * line_width + (zero_based % line_bases)
                handle.seek(byte_position)
                observed.append(handle.read(1).decode("ascii").upper())
            observed_reference = "".join(observed)
            verified = observed_reference == reference
            item.update(
                reference_allele_verified=verified,
                observed_reference=observed_reference,
                reference_validation_status="verified" if verified else "mismatch",
            )
            if not verified:
                blockers.append(f"reference_mismatch:{item['variant']}")
                result["omitted"].append({**item, "omission_reason": "reference_allele_mismatch"})
            else:
                validated_selected.append(item)
    result["selected"] = validated_selected
    result["blockers"] = sorted(set(blockers))
    result["status"] = "available" if result["selected"] and not blockers else "blocked"
    return result


def prepare_alphagenome_disclosure(
    selection: dict[str, Any],
    *,
    ontology_terms: Iterable[str],
    modalities: Iterable[str],
    sequence_length: int,
) -> dict[str, Any]:
    terms = list(dict.fromkeys(_text(term) for term in ontology_terms if _text(term)))
    requested = list(dict.fromkeys(_text(item).upper() for item in modalities if _text(item)))
    blockers = list(selection.get("blockers") or [])
    if not 1 <= len(terms) <= MAX_ONTOLOGY_TERMS:
        blockers.append("ontology_terms_must_contain_1_to_5_items")
    if not requested or any(item not in SUPPORTED_ALPHAGENOME_MODALITIES for item in requested):
        blockers.append("unsupported_or_missing_modalities")
    if sequence_length not in SUPPORTED_SEQUENCE_LENGTHS:
        blockers.append("unsupported_sequence_length")
    variants = []
    if not blockers:
        for item in selection.get("selected", []):
            center = int(item["position"]) - 1
            start = max(0, center - sequence_length // 2)
            variants.append(
                {
                    "variant": item["variant"],
                    "assembly": "GRCh38",
                    "chromosome": item["chromosome"],
                    "position_1_based": item["position"],
                    "reference": item["reference"],
                    "alternate": item["alternate"],
                    "model_interval_0_based_half_open": {
                        "start": start,
                        "end": start + sequence_length,
                        "length": sequence_length,
                    },
                    "reference_allele_verified": bool(item.get("reference_allele_verified", False)),
                }
            )
    payload = {
        "model_id": "alphagenome-api",
        "assembly": "GRCh38",
        "variants": variants,
        "ontology_terms": terms,
        "ontology_handling": "worker_local_result_filter_not_provider_request_parameter",
        "modalities": requested,
        "sequence_length": sequence_length,
        "external_transfer": True,
    }
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
    return {
        "status": "ready" if not blockers and variants else "blocked",
        "blockers": sorted(set(blockers or (["no_selected_variants"] if not variants else []))),
        "payload": payload,
        "payload_sha256": digest,
        "selected_variant_count": len(variants),
        "omitted_variant_count": len(selection.get("omitted", [])),
    }
