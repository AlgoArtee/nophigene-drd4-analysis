"""Coverage for the compressed gene-data archive format."""

from __future__ import annotations

import json
import zipfile
from pathlib import Path

import pandas as pd

import src.analysis as analysis
from src.helper_functions.filter_manifest_region import sanitize_gene_name_for_filename


def _write_zip_member(bundle: zipfile.ZipFile, name: str, payload: object | str) -> None:
    if isinstance(payload, str):
        bundle.writestr(name, payload)
    else:
        bundle.writestr(name, json.dumps(payload))


def _install_test_bundle(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(analysis, "GENE_DATA_DIR", tmp_path)
    monkeypatch.setattr(analysis, "GENE_DATA_BUNDLE_PATH", tmp_path / "gene_data_bundle.zip")
    analysis.clear_gene_data_bundle_cache()

    interpretation = {
        "database_name": "NophiGene TEST1 Interpretation Database",
        "gene_context": {
            "gene_name": "TEST1",
            "gene_region": {"start": 101, "end": 202},
            "variant_effect_overview": ["archive interpretation"],
            "relevant_methylation_probe_ids": ["cg00000001"],
        },
        "variant_records": [],
    }
    population = {
        "database_name": "NophiGene TEST1 Population Database",
        "gene_population_patterns": [{"summary": "archive population"}],
    }
    manifest = "IlmnID,CHR,MAPINFO,UCSC_RefGene_Name\ncg00000001,1,101,TEST1\n"

    with zipfile.ZipFile(analysis.GENE_DATA_BUNDLE_PATH, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
        _write_zip_member(bundle, "test1_interpretation_db.json", interpretation)
        _write_zip_member(bundle, "test1_population_db.json", population)
        _write_zip_member(bundle, "TEST1_epigenetics_hg19.csv", manifest)


def test_gene_data_loaders_fall_back_to_archive(monkeypatch, tmp_path: Path) -> None:
    _install_test_bundle(monkeypatch, tmp_path)

    knowledge_base = analysis.load_gene_interpretation_database("TEST1")
    population_database = analysis.load_gene_population_database("test1")
    manifest = analysis.load_gene_epigenetics_manifest("TEST1")

    assert knowledge_base is not None
    assert population_database is not None
    assert manifest is not None
    assert knowledge_base["gene_context"]["gene_region"]["start"] == 101
    assert population_database["gene_population_patterns"][0]["summary"] == "archive population"
    assert manifest["IlmnID"].tolist() == ["cg00000001"]


def test_loose_gene_data_overrides_archive(monkeypatch, tmp_path: Path) -> None:
    _install_test_bundle(monkeypatch, tmp_path)

    loose_interpretation = {
        "database_name": "NophiGene TEST1 Interpretation Database",
        "gene_context": {
            "gene_name": "TEST1",
            "gene_region": {"start": 303, "end": 404},
            "variant_effect_overview": ["loose interpretation"],
            "relevant_methylation_probe_ids": ["cg00000002"],
        },
        "variant_records": [],
    }
    (tmp_path / "test1_interpretation_db.json").write_text(json.dumps(loose_interpretation), encoding="utf-8")
    (tmp_path / "TEST1_epigenetics_hg19.csv").write_text(
        "IlmnID,CHR,MAPINFO,UCSC_RefGene_Name\ncg00000002,1,303,TEST1\n",
        encoding="utf-8",
    )

    knowledge_base = analysis.load_gene_interpretation_database("TEST1")
    manifest = analysis.load_gene_epigenetics_manifest("TEST1")

    assert knowledge_base is not None
    assert manifest is not None
    assert knowledge_base["gene_context"]["gene_region"]["start"] == 303
    assert manifest["IlmnID"].tolist() == ["cg00000002"]


def test_missing_gene_data_still_returns_none(monkeypatch, tmp_path: Path) -> None:
    _install_test_bundle(monkeypatch, tmp_path)

    assert analysis.load_gene_interpretation_database("NOPE") is None
    assert analysis.load_gene_population_database("NOPE") is None
    assert analysis.load_gene_epigenetics_manifest("NOPE") is None


def test_gene_data_bundle_contains_valid_generated_triplets() -> None:
    required_csv_columns = {"IlmnID", "CHR", "MAPINFO", "UCSC_RefGene_Name"}
    bundle_path = analysis.GENE_DATA_BUNDLE_PATH

    assert bundle_path.exists()
    with zipfile.ZipFile(bundle_path) as bundle:
        members = {info.filename for info in bundle.infolist() if not info.is_dir()}
        interpretation_members = sorted(name for name in members if name.endswith("_interpretation_db.json"))
        assert interpretation_members

        for member in interpretation_members:
            payload = json.loads(bundle.read(member).decode("utf-8"))
            gene_name = payload["gene_context"]["gene_name"]
            file_gene_name = sanitize_gene_name_for_filename(gene_name)
            assert f"{file_gene_name.lower()}_population_db.json" in members
            csv_name = f"{file_gene_name}_epigenetics_hg19.csv"
            assert csv_name in members

            csv_header = pd.read_csv(bundle.open(csv_name), nrows=0).columns
            assert required_csv_columns.issubset(set(csv_header))

        for member in sorted(name for name in members if name.endswith(".json")):
            json.loads(bundle.read(member).decode("utf-8"))
