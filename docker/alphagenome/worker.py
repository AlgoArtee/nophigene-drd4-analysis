"""Minimal signed-job AlphaGenome worker with a fixed provider endpoint."""

from __future__ import annotations

import hashlib
import hmac
import json
import math
import os
import time
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from typing import Any

import pyzipper

from prepare_request import prepare_request


WORK_ROOT = Path(os.environ.get("NOPHIGENE_MODEL_WORK_ROOT", "/work"))
CREDENTIAL_ROOT = Path(os.environ.get("NOPHIGENE_MODEL_CREDENTIAL_ROOT", "/credentials"))
QUEUE = WORK_ROOT / "queue"
JOBS = WORK_ROOT / "jobs"
PROVIDER_ADDRESS = "dns:///gdmscience.googleapis.com:443"
CONTRACT_VERSION = "1.0"
REFERENCE_FASTA = Path(os.environ.get("NOPHIGENE_HG38_REFERENCE_FASTA", "/reference/hg38.analysisSet.fa"))
WORKER_CONTAINER_VERSION = os.environ.get(
    "NOPHIGENE_MODEL_WORKER_IMAGE", "nophigene/alphagenome-runner:0.8.0"
)


def now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def canonical_json(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, default=str).encode("utf-8")


def secret(path_env: str) -> bytes:
    path = os.environ.get(path_env, "")
    if not path:
        raise RuntimeError(f"Required secret file is missing: {path_env}")
    value = Path(path).read_bytes().strip()
    if len(value) < 32:
        raise RuntimeError(f"Required secret is invalid: {path_env}")
    return value


def verify(envelope: dict[str, Any]) -> dict[str, Any]:
    payload = envelope.get("payload")
    signature = str(envelope.get("signature") or "")
    if not isinstance(payload, dict) or envelope.get("algorithm") != "HMAC-SHA256":
        raise RuntimeError("Invalid signed model-job envelope.")
    expected = hmac.new(secret("NOPHIGENE_MODEL_RUNNER_KEY_FILE"), canonical_json(payload), hashlib.sha256).hexdigest()
    if not hmac.compare_digest(signature, expected):
        raise RuntimeError("Model-job envelope failed authentication.")
    return payload


def sign(payload: dict[str, Any]) -> dict[str, Any]:
    signature = hmac.new(secret("NOPHIGENE_MODEL_RUNNER_KEY_FILE"), canonical_json(payload), hashlib.sha256).hexdigest()
    return {"payload": payload, "signature": signature, "algorithm": "HMAC-SHA256"}


def api_key() -> str:
    archive_path = CREDENTIAL_ROOT / "alphagenome-api.aes.zip"
    if not archive_path.is_file():
        raise RuntimeError("AlphaGenome credential is not configured.")
    with pyzipper.AESZipFile(archive_path, "r") as archive:
        archive.setpassword(secret("NOPHIGENE_MODEL_CREDENTIAL_KEY_FILE"))
        value = archive.read("credential.txt").decode("utf-8").strip()
    if not value:
        raise RuntimeError("AlphaGenome credential is empty.")
    return value


def redacted_error(error: Exception) -> str:
    """Return a bounded error string with the configured credential removed."""
    message = str(error)
    try:
        configured_key = api_key()
    except Exception:
        configured_key = ""
    if configured_key:
        message = message.replace(configured_key, "[redacted-provider-credential]")
    return message[:500]


def write_status(job_id: str, **updates: Any) -> dict[str, Any]:
    path = JOBS / job_id / "status.json"
    state = json.loads(path.read_text(encoding="utf-8"))
    state.update(updates)
    temporary = path.with_suffix(".tmp")
    temporary.write_bytes(canonical_json(state))
    temporary.replace(path)
    return state


def output_type_name(value: Any) -> str:
    return str(value).split(".")[-1].upper()


def scorer_matches(scorer: Any, modalities: set[str]) -> bool:
    try:
        serialized = f"{type(scorer).__name__} {scorer.to_proto()}".upper()
    except Exception:
        serialized = repr(scorer).upper()
    normalized = "".join(character for character in serialized if character.isalnum())
    return any("".join(character for character in modality if character.isalnum()) in normalized for modality in modalities)


def reference_bases(chromosome: str, position_1_based: int, length: int) -> str:
    """Read reference bases directly from the independently mounted FASTA."""
    index_path = Path(f"{REFERENCE_FASTA}.fai")
    if not REFERENCE_FASTA.is_file() or not index_path.is_file():
        raise RuntimeError("The worker GRCh38 reference FASTA or index is unavailable.")
    requested = str(chromosome).removeprefix("chr").upper()
    metadata = None
    for line in index_path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) < 5 or fields[0].removeprefix("chr").upper() != requested:
            continue
        metadata = tuple(int(value) for value in fields[1:5])
        break
    if metadata is None:
        raise RuntimeError(f"Reference chromosome is unavailable: {chromosome}")
    sequence_length, offset, line_bases, line_width = metadata
    bases: list[str] = []
    with REFERENCE_FASTA.open("rb") as handle:
        for delta in range(length):
            zero_based = position_1_based - 1 + delta
            if zero_based < 0 or zero_based >= sequence_length:
                break
            byte_position = offset + (zero_based // line_bases) * line_width + (zero_based % line_bases)
            handle.seek(byte_position)
            bases.append(handle.read(1).decode("ascii").upper())
    return "".join(bases)


def tidy_rows(frame: Any, *, requested_modalities: set[str], ontology_terms: set[str]) -> list[dict[str, Any]]:
    if frame is None:
        return []
    rows: list[dict[str, Any]] = []
    for raw in frame.to_dict(orient="records"):
        output_type = output_type_name(raw.get("output_type"))
        ontology = str(raw.get("ontology_curie") or "")
        if output_type and output_type not in requested_modalities:
            continue
        if ontology_terms and ontology and ontology not in ontology_terms:
            continue
        variant = str(raw.get("variant_id") or raw.get("variant") or "")
        scorer = str(raw.get("variant_scorer") or raw.get("interval_scorer") or "")
        track = str(raw.get("track_name") or "")
        gene_context = str(raw.get("gene_id") or raw.get("gene_name") or "")
        entity_key = variant
        output_identity = "|".join((output_type, scorer, ontology, track, gene_context))
        output_name = f"{output_type}:{hashlib.sha256(output_identity.encode()).hexdigest()[:16]}"
        raw_score = raw.get("raw_score")
        try:
            raw_score = float(raw_score)
        except (TypeError, ValueError):
            raw_score = None
        if raw_score is not None and not math.isfinite(raw_score):
            raw_score = None
        quantile = raw.get("quantile_score")
        try:
            quantile = float(quantile)
        except (TypeError, ValueError):
            quantile = None
        if quantile is not None and not math.isfinite(quantile):
            quantile = None
        rows.append(
            {
                "entity_type": "observed_variant",
                "entity_key": entity_key,
                "output_name": output_name,
                "variant": variant,
                "output_type": output_type,
                "scorer": scorer,
                "track_name": track,
                "gene_id": raw.get("gene_id"),
                "gene_name": raw.get("gene_name"),
                "gene_type": raw.get("gene_type"),
                "gene_strand": raw.get("gene_strand"),
                "ontology_curie": ontology,
                "gtex_tissue": raw.get("gtex_tissue"),
                "biosample_name": raw.get("biosample_name"),
                "biosample_type": raw.get("biosample_type"),
                "transcription_factor": raw.get("transcription_factor"),
                "histone_mark": raw.get("histone_mark"),
                "raw_score": raw_score,
                "quantile_score": quantile,
                "validated_label": "molecular_variant_effect_score",
                "applicability": {
                    "ontology_curie": ontology,
                    "biosample_name": raw.get("biosample_name"),
                    "gene_id": raw.get("gene_id"),
                    "gene_name": raw.get("gene_name"),
                    "track_name": track,
                    "output_type": output_type,
                    "scorer": scorer,
                    "quantile_score": quantile,
                },
                "limitations": [
                    "AlphaGenome molecular prediction; not a pathogenicity probability, statistical significance measurement, diagnosis, or disease-risk estimate.",
                    "Score sign and magnitude require scorer-specific interpretation.",
                ],
            }
        )
    return rows


def sanitized_metadata(model: Any) -> dict[str, Any]:
    from alphagenome.models import dna_client

    metadata = model.output_metadata(dna_client.Organism.HOMO_SAPIENS).concatenate()
    rows = []
    for item in metadata.to_dict(orient="records"):
        curie = str(item.get("ontology_curie") or "")
        if not curie:
            continue
        rows.append(
            {
                "ontology_curie": curie,
                "biosample_name": str(item.get("biosample_name") or ""),
                "biosample_type": str(item.get("biosample_type") or ""),
                "output_type": output_type_name(item.get("output_type")),
            }
        )
    unique = {json.dumps(item, sort_keys=True): item for item in rows}
    return {
        "sdk_version": version("alphagenome"),
        "sdk_source_commit": "71a6beb8c30832f121309a81c2530efa5af7986a",
        "worker_container_version": WORKER_CONTAINER_VERSION,
        "provider_model_revision": None,
        "retrieved_at": now(),
        "ontology_terms": sorted(unique.values(), key=lambda item: (item["biosample_name"], item["ontology_curie"], item["output_type"])),
    }


def execute(payload: dict[str, Any]) -> dict[str, Any]:
    if payload.get("contract_version") != CONTRACT_VERSION or payload.get("model_id") != "alphagenome-api":
        raise RuntimeError("Unsupported model runner contract or model identifier.")
    inputs = payload.get("inputs")
    if not isinstance(inputs, dict):
        raise RuntimeError("Model inputs are invalid.")
    if (JOBS / str(payload.get("job_id") or "") / "cancel.requested").exists():
        raise RuntimeError("Model job cancelled before provider access.")
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers

    try:
        client = dna_client.create(api_key(), address=PROVIDER_ADDRESS, timeout=30)
    except Exception as exc:
        raise RuntimeError(f"AlphaGenome client initialization failed: {redacted_error(exc)}") from None
    if inputs.get("operation") == "metadata":
        return {"kind": "metadata", **sanitized_metadata(client), "predictions": [], "failures": []}

    modalities = {str(item).upper() for item in inputs.get("requested_modalities", [])}
    ontology_terms = {str(item) for item in inputs.get("ontology_terms", [])}
    recommended = variant_scorers.get_recommended_scorers(dna_client.Organism.HOMO_SAPIENS.to_proto())
    selected_scorers = [scorer for scorer in recommended if scorer_matches(scorer, modalities)]
    if not selected_scorers:
        raise RuntimeError("No recommended AlphaGenome scorers match the requested modalities.")
    predictions: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    variants = inputs.get("variants") or []
    for index, item in enumerate(variants, start=1):
        if (JOBS / payload["job_id"] / "cancel.requested").exists():
            raise RuntimeError("Model job cancelled.")
        if item.get("assembly") != "GRCh38" or item.get("reference_allele_verified") is not True:
            failures.append({"variant": item.get("variant"), "code": "reference_not_verified"})
            continue
        try:
            canonical = prepare_request(
                variant=str(item.get("variant") or ""),
                assembly=str(item.get("assembly") or ""),
                modalities=modalities,
                ontology_terms=ontology_terms,
                sequence_length=int(inputs.get("sequence_length") or 0),
                external_transfer_acknowledged=bool(inputs.get("explicit_transfer_consent")),
            )
            interval_data = canonical["model_interval_0_based_half_open"]
            if interval_data != item.get("model_interval_0_based_half_open"):
                raise RuntimeError("The signed interval does not match the canonical request.")
            observed_reference = reference_bases(
                str(item["chromosome"]), int(item["position_1_based"]), len(str(item["reference"]))
            )
            if observed_reference != str(item["reference"]).upper():
                failures.append(
                    {
                        "variant": item.get("variant"),
                        "code": "worker_reference_mismatch",
                        "observed_reference": observed_reference,
                    }
                )
                continue
            interval = genome.Interval(item["chromosome"], int(interval_data["start"]), int(interval_data["end"]))
            variant = genome.Variant(
                chromosome=item["chromosome"],
                position=int(item["position_1_based"]),
                reference_bases=item["reference"],
                alternate_bases=item["alternate"],
            )
            scores = client.score_variant(interval=interval, variant=variant, variant_scorers=selected_scorers)
            frame = variant_scorers.tidy_scores(scores)
            variant_rows = tidy_rows(frame, requested_modalities=modalities, ontology_terms=ontology_terms)
            for row in variant_rows:
                row["variant"] = item["variant"]
                row["entity_key"] = item["variant"]
                row["model_interval_0_based_half_open"] = interval_data
            predictions.extend(variant_rows)
            write_status(payload["job_id"], progress_percent=max(5, int(index / len(variants) * 95)))
        except Exception as exc:
            failures.append({"variant": item.get("variant"), "code": "provider_variant_failed", "message": redacted_error(exc)})
    return {
        "kind": "predictions",
        "model_id": "alphagenome-api",
        "sdk_version": version("alphagenome"),
        "sdk_source_commit": "71a6beb8c30832f121309a81c2530efa5af7986a",
        "worker_container_version": WORKER_CONTAINER_VERSION,
        "provider_model_revision": None,
        "provider_endpoint": PROVIDER_ADDRESS,
        "requested_modalities": sorted(modalities),
        "ontology_terms": sorted(ontology_terms),
        "retrieved_at": now(),
        "external_transfer_occurred": True,
        "predictions": predictions,
        "failures": failures,
        "empty_output": not predictions and not failures,
        "limitations": [
            "Raw and quantile scores are molecular model outputs, not statistical significance measurements or clinical classifications.",
            "The public API did not expose an immutable provider model revision for this request.",
        ],
    }


def recover_completed_result(job_id: str) -> bool:
    """Finalize a signed result left behind by a worker restart."""
    result_path = JOBS / job_id / "result.json"
    if not result_path.is_file():
        return False
    envelope = json.loads(result_path.read_text(encoding="utf-8"))
    result = verify(envelope)
    claimed_checksum = str(result.get("result_checksum_sha256") or "")
    checksum_payload = {key: value for key, value in result.items() if key != "result_checksum_sha256"}
    actual_checksum = hashlib.sha256(canonical_json(checksum_payload)).hexdigest()
    if not claimed_checksum or not hmac.compare_digest(claimed_checksum, actual_checksum):
        raise RuntimeError("Recovered result checksum validation failed.")
    prediction_count = len(result.get("predictions") or [])
    failure_count = len(result.get("failures") or [])
    final_status = "partial" if prediction_count and failure_count else "succeeded"
    if not prediction_count and failure_count:
        final_status = "failed"
    write_status(
        job_id,
        status=final_status,
        stage="completed_after_worker_restart" if final_status in {"succeeded", "partial"} else "failed",
        progress_percent=100,
        finished_at=now(),
        result_checksum_sha256=claimed_checksum,
        prediction_count=prediction_count,
        failure_count=failure_count,
        empty_output=bool(result.get("empty_output")),
        error=None if final_status != "failed" else {"code": "all_variants_failed", "message": "Every selected variant failed."},
    )
    return True


def process(path: Path) -> None:
    envelope = json.loads(path.read_text(encoding="utf-8"))
    payload = verify(envelope)
    job_id = str(payload.get("job_id") or "")
    if len(job_id) != 32 or (JOBS / job_id).resolve().parent != JOBS.resolve():
        raise RuntimeError("Unsafe model job identifier.")
    if recover_completed_result(job_id):
        path.unlink(missing_ok=True)
        return
    write_status(job_id, status="running", stage="provider_inference", started_at=now(), progress_percent=5)
    try:
        result = execute(payload)
        result_checksum = hashlib.sha256(canonical_json(result)).hexdigest()
        result_path = JOBS / job_id / "result.json"
        temporary_result_path = result_path.with_suffix(".tmp")
        temporary_result_path.write_bytes(canonical_json(sign({**result, "result_checksum_sha256": result_checksum})))
        temporary_result_path.replace(result_path)
        prediction_count = len(result.get("predictions", []))
        failure_count = len(result.get("failures", []))
        empty_output = bool(result.get("empty_output"))
        final_status = "partial" if prediction_count and failure_count else "succeeded"
        if not prediction_count and failure_count:
            final_status = "failed"
        write_status(
            job_id,
            status=final_status,
            stage="completed" if final_status in {"succeeded", "partial"} else "failed",
            progress_percent=100,
            finished_at=now(),
            result_checksum_sha256=result_checksum,
            prediction_count=prediction_count,
            failure_count=failure_count,
            empty_output=empty_output,
            error=None if final_status != "failed" else {"code": "all_variants_failed", "message": "Every selected variant failed."},
        )
    except Exception as exc:
        write_status(
            job_id,
            status="cancelled" if "cancelled" in str(exc).casefold() else "failed",
            stage="failed",
            progress_percent=100,
            finished_at=now(),
            error={"code": "model_execution_failed", "message": redacted_error(exc)},
        )
    path.unlink(missing_ok=True)


def heartbeat(status: str = "ready", reason: str = "") -> None:
    payload = {
        "status": status,
        "reason": reason,
        "updated_at": now(),
        "worker": "alphagenome",
        "sdk_version": version("alphagenome"),
        "worker_container_version": WORKER_CONTAINER_VERSION,
        "provider_endpoint": PROVIDER_ADDRESS,
    }
    (WORK_ROOT / "worker-heartbeat.json").write_bytes(canonical_json(payload))


def main() -> None:
    QUEUE.mkdir(parents=True, exist_ok=True)
    JOBS.mkdir(parents=True, exist_ok=True)
    while True:
        heartbeat()
        queued = sorted(QUEUE.glob("*.json"))
        if not queued:
            time.sleep(2)
            continue
        try:
            process(queued[0])
        except Exception as exc:
            heartbeat("degraded", redacted_error(exc)[:300])
            time.sleep(2)


if __name__ == "__main__":
    main()
