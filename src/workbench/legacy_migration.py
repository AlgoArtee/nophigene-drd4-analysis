"""Dry-run-first import and encrypted retirement of legacy stores."""

from __future__ import annotations

import hashlib
import json
import zipfile
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

from sqlalchemy import select
from sqlalchemy.orm import Session

from .artifacts import ArtifactStore
from .audit import append_audit_event
from .models import EvidenceRecord, EvidenceSnapshot, Run
from .reporting import build_canonical_report


@dataclass
class LegacyCandidate:
    path: str
    kind: str
    status: str
    checksum_sha256: str = ""
    gene: str = ""
    schema_version: str = ""
    excluded_fields: list[str] = field(default_factory=list)
    error: str = ""


@dataclass
class LegacyMigrationReport:
    candidates: list[LegacyCandidate] = field(default_factory=list)

    @property
    def counts(self) -> dict[str, int]:
        result: dict[str, int] = {}
        for candidate in self.candidates:
            result[candidate.status] = result.get(candidate.status, 0) + 1
        return result

    def to_dict(self) -> dict[str, Any]:
        return {"counts": self.counts, "candidates": [asdict(item) for item in self.candidates]}


def _digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def _gene_from_filename(path: Path, suffix: str) -> str:
    return path.name.removesuffix(suffix).replace("_", "-").upper()


def _legacy_reference_links(value: Any) -> list[dict[str, str]]:
    """Extract only explicit bibliographic/linkout objects from a legacy payload."""
    found: dict[tuple[str, str], dict[str, str]] = {}

    def visit(item: Any) -> None:
        if isinstance(item, dict):
            url = str(item.get("url") or "").strip()
            label = str(item.get("label") or item.get("title") or "").strip()
            if url.startswith(("https://", "http://")) and label:
                found[(url.casefold(), label.casefold())] = {"url": url, "label": label}
            for child in item.values():
                visit(child)
        elif isinstance(item, list):
            for child in item:
                visit(child)

    visit(value)
    return list(found.values())


def _legacy_source_key(url: str) -> str:
    host = (urlparse(url).hostname or "legacy-linkout").lower()
    aliases = {
        "pubmed.ncbi.nlm.nih.gov": "pubmed",
        "www.ncbi.nlm.nih.gov": "ncbi",
        "www.uniprot.org": "uniprot",
        "grch37.ensembl.org": "ensembl",
        "www.ensembl.org": "ensembl",
    }
    return aliases.get(host, host.removeprefix("www.")[:64])


def inventory_legacy_stores(project_root: Path) -> LegacyMigrationReport:
    root = Path(project_root)
    report = LegacyMigrationReport()
    for path in sorted((root / "results").rglob("*.json")) if (root / "results").exists() else []:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            if not isinstance(payload, dict):
                raise ValueError("root is not an object")
            excluded = [key for key in ("predictive_theses", "synthesis", "phenotype_prediction") if key in payload]
            gene = str(payload.get("gene") or payload.get("gene_name") or payload.get("run", {}).get("gene") or "").upper()
            report.candidates.append(
                LegacyCandidate(
                    path=str(path), kind="report_json", status="parseable", checksum_sha256=_digest(path),
                    gene=gene, schema_version=str(payload.get("schema_version") or "legacy"), excluded_fields=excluded,
                )
            )
        except Exception as exc:
            report.candidates.append(LegacyCandidate(path=str(path), kind="report_json", status="quarantined", error=str(exc)))
    central_csv = root / "results" / "general_gene_analysis_database.csv"
    if central_csv.exists():
        report.candidates.append(
            LegacyCandidate(
                path=str(central_csv), kind="denormalized_central_csv", status="archive_only",
                checksum_sha256=_digest(central_csv), excluded_fields=["duplicated_gene_level_methylation_context"],
            )
        )
    active_gene_data = root / "src" / "gene_data"
    archived_loose_gene_data = root / "version1" / "gene_data" / "loose"
    gene_data_roots = [path for path in (active_gene_data, archived_loose_gene_data) if path.exists()]
    for gene_data in gene_data_roots:
        for path in sorted(gene_data.glob("*_synthesis.json")):
            report.candidates.append(
                LegacyCandidate(
                    path=str(path), kind="synthetic_synthesis", status="archive_only",
                    checksum_sha256=_digest(path), gene=_gene_from_filename(path, "_synthesis.json"),
                    excluded_fields=["entire_file"],
                )
            )
    for suffix, kind in (
        ("_interpretation_db.json", "curated_interpretation_references"),
        ("_population_db.json", "curated_population_references"),
    ):
        for gene_data in gene_data_roots:
            for path in sorted(gene_data.glob(f"*{suffix}")):
                try:
                    payload = json.loads(path.read_text(encoding="utf-8"))
                    if not isinstance(payload, dict):
                        raise ValueError("root is not an object")
                    report.candidates.append(
                        LegacyCandidate(
                            path=str(path), kind=kind, status="evidence_import_candidate",
                            checksum_sha256=_digest(path), gene=_gene_from_filename(path, suffix),
                            schema_version=str(payload.get("version") or "legacy"),
                            excluded_fields=["concrete_variant_prediction", "predictive_theses", "synthetic_cases"],
                        )
                    )
                except Exception as exc:
                    report.candidates.append(
                        LegacyCandidate(path=str(path), kind=kind, status="quarantined", error=str(exc))
                    )
    for path in sorted(active_gene_data.glob("*.zip")) if active_gene_data.exists() else []:
        report.candidates.append(
            LegacyCandidate(
                path=str(path), kind="gene_bundle", status="evidence_import_candidate",
                checksum_sha256=_digest(path), excluded_fields=["*_synthesis.json", "synthetic_cases"],
            )
        )
    return report


def import_legacy_reports(
    session: Session,
    migration_report: LegacyMigrationReport,
    *,
    artifact_store: ArtifactStore,
) -> dict[str, Any]:
    imported = 0
    skipped = 0
    evidence_references_imported = 0
    quarantined: list[dict[str, str]] = []
    for candidate in migration_report.candidates:
        if candidate.kind in {"curated_interpretation_references", "curated_population_references"} and candidate.status == "evidence_import_candidate":
            path = Path(candidate.path)
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
                links = _legacy_reference_links(payload)
                if not links:
                    skipped += 1
                    continue
                existing = session.scalar(
                    select(EvidenceSnapshot).where(EvidenceSnapshot.checksum_sha256 == candidate.checksum_sha256)
                )
                if existing:
                    skipped += 1
                    continue
                raw_artifact = artifact_store.put_bytes(
                    session,
                    path.read_bytes(),
                    kind="legacy_curated_reference_payload",
                    filename=path.name,
                    media_type="application/json",
                    sensitive=False,
                )
                snapshot_id = hashlib.sha256(("legacy-evidence:" + candidate.checksum_sha256).encode("utf-8")).hexdigest()[:32]
                session.add(
                    EvidenceSnapshot(
                        id=snapshot_id,
                        checksum_sha256=candidate.checksum_sha256,
                        refresh_policy="immutable_legacy_import",
                        query_context={
                            "gene": candidate.gene,
                            "legacy_path": candidate.path,
                            "import_state": "linkouts_not_assessed",
                        },
                    )
                )
                session.flush()
                for link in links:
                    record_id = hashlib.sha256((link["url"] + "\n" + link["label"]).encode("utf-8")).hexdigest()
                    session.add(
                        EvidenceRecord(
                            snapshot_id=snapshot_id,
                            source_key=_legacy_source_key(link["url"]),
                            source_release=candidate.schema_version,
                            record_id=record_id,
                            evidence_type="bibliographic_linkout",
                            entity_type="gene",
                            entity_key=candidate.gene,
                            assertion=link["label"],
                            evidence_level="unverified_legacy_reference",
                            context={"url": link["url"], "checked": False, "legacy_kind": candidate.kind},
                            citations=[{"url": link["url"], "label": link["label"]}],
                            license_status="not_assessed",
                            raw_checksum_sha256=candidate.checksum_sha256,
                            raw_artifact_id=raw_artifact.id,
                            status="not_assessed",
                        )
                    )
                    evidence_references_imported += 1
                append_audit_event(
                    session,
                    "legacy_evidence_references_imported",
                    entity_type="evidence_snapshot",
                    entity_id=snapshot_id,
                    payload={
                        "source_checksum_sha256": candidate.checksum_sha256,
                        "reference_count": len(links),
                        "status": "not_assessed",
                        "excluded_fields": candidate.excluded_fields,
                    },
                )
                imported += 1
            except Exception as exc:
                quarantined.append({"path": candidate.path, "error": str(exc)})
            continue
        if candidate.kind != "report_json" or candidate.status != "parseable":
            continue
        path = Path(candidate.path)
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload.pop("predictive_theses", None)
            payload.pop("synthesis", None)
            canonical = payload if payload.get("schema_version") == "3.0" else build_canonical_report(payload)
            run_id = hashlib.sha256((candidate.checksum_sha256 + candidate.path).encode("utf-8")).hexdigest()[:32]
            if session.get(Run, run_id):
                skipped += 1
                continue
            run = Run(
                id=run_id,
                schema_version="3.0",
                status="succeeded",
                stage="legacy_import",
                progress_percent=100,
                genes=[canonical.get("run", {}).get("gene")] if canonical.get("run", {}).get("gene") else [],
                context=canonical.get("run", {}).get("sample_context", {}),
                configuration={"legacy_path": candidate.path, "legacy_schema": candidate.schema_version},
            )
            session.add(run)
            session.flush()
            artifact_store.put_bytes(
                session,
                json.dumps(canonical, ensure_ascii=False, indent=2).encode("utf-8"),
                kind="canonical_legacy_report",
                run_id=run_id,
                filename="report.json",
                media_type="application/json",
                sensitive=True,
            )
            append_audit_event(
                session,
                "legacy_report_imported",
                entity_type="run",
                entity_id=run_id,
                payload={"source_checksum_sha256": candidate.checksum_sha256, "excluded_fields": candidate.excluded_fields},
            )
            imported += 1
        except Exception as exc:
            quarantined.append({"path": candidate.path, "error": str(exc)})
    return {
        "imported": imported,
        "evidence_references_imported": evidence_references_imported,
        "skipped": skipped,
        "quarantined": quarantined,
    }


def create_encrypted_legacy_archive(
    migration_report: LegacyMigrationReport,
    output_path: Path,
    *,
    password: str,
) -> dict[str, Any]:
    if not password:
        raise ValueError("An archive password is required.")
    try:
        import pyzipper
    except ImportError as exc:
        raise RuntimeError("Encrypted legacy archives require pyzipper.") from exc
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        raise FileExistsError(f"Refusing to overwrite immutable legacy archive: {output_path}")
    manifest = migration_report.to_dict()
    with pyzipper.AESZipFile(
        output_path, "w", compression=zipfile.ZIP_DEFLATED, encryption=pyzipper.WZ_AES
    ) as archive:
        archive.setpassword(password.encode("utf-8"))
        archive.writestr("migration-manifest.json", json.dumps(manifest, indent=2).encode("utf-8"))
        for candidate in migration_report.candidates:
            path = Path(candidate.path)
            if path.is_file():
                archive.write(path, f"legacy/{candidate.kind}/{candidate.checksum_sha256[:12]}-{path.name}")
    return {"path": str(output_path), "checksum_sha256": _digest(output_path), "file_count": len(migration_report.candidates)}
