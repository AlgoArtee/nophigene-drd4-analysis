"""Coverage for sharded curated-lite protein-coding gene data."""

from __future__ import annotations

import gzip
import json
import zipfile
from pathlib import Path

import pandas as pd

import src.analysis as analysis
from scripts.generate_bulk_protein_coding_gene_data import build_bulk_gene_data


def _write_zip_member(bundle: zipfile.ZipFile, name: str, payload: object | str) -> None:
    if isinstance(payload, str):
        bundle.writestr(name, payload)
    else:
        bundle.writestr(name, json.dumps(payload))


def _minimal_interpretation(gene_name: str, start: int) -> dict[str, object]:
    return {
        "database_name": f"NophiGene {gene_name} Interpretation Database",
        "gene_context": {
            "gene_name": gene_name,
            "chromosome": "1",
            "gene_region": {"start": start, "end": start + 10},
            "promoter_review_region": {"start": max(1, start - 1000), "end": start - 1},
            "recommended_promoter_plus_gene_region": f"1:{max(1, start - 1000)}-{start + 10}",
            "variant_effect_overview": ["test"],
            "relevant_methylation_probe_ids": [],
        },
        "variant_records": [],
    }


def _install_bulk_test_data(monkeypatch, tmp_path: Path) -> None:
    gene_data_dir = tmp_path / "gene_data"
    shard_dir = gene_data_dir / "bulk_gene_data_shards"
    shard_dir.mkdir(parents=True)
    monkeypatch.setattr(analysis, "GENE_DATA_DIR", gene_data_dir)
    monkeypatch.setattr(analysis, "GENE_DATA_BUNDLE_PATH", gene_data_dir / "gene_data_bundle.zip")
    monkeypatch.setattr(analysis, "GENE_DATA_INDEX_PATH", gene_data_dir / "gene_data_index.json")
    analysis.clear_gene_data_bundle_cache()

    with zipfile.ZipFile(analysis.GENE_DATA_BUNDLE_PATH, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
        _write_zip_member(bundle, "curated1_interpretation_db.json", _minimal_interpretation("CURATED1", 11))
        _write_zip_member(bundle, "curated1_population_db.json", {"database_name": "curated population"})
        _write_zip_member(bundle, "CURATED1_epigenetics_hg19.csv", "IlmnID,CHR,MAPINFO,UCSC_RefGene_Name\ncg1,1,11,CURATED1\n")

    shard_name = "gene_data_bulk_shard_00.zip"
    with zipfile.ZipFile(shard_dir / shard_name, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
        _write_zip_member(bundle, "bulk1_interpretation_db.json", _minimal_interpretation("BULK1", 101))
        _write_zip_member(bundle, "bulk1_population_db.json", {"database_name": "bulk population", "gene_population_patterns": []})
        _write_zip_member(bundle, "BULK1_epigenetics_hg19.csv", "IlmnID,CHR,MAPINFO,UCSC_RefGene_Name\n")
        _write_zip_member(bundle, "curated1_interpretation_db.json", _minimal_interpretation("CURATED1", 999))

    index = {
        "format_version": "test",
        "genes": {
            "BULK1": {
                "shard": shard_name,
                "files": {
                    "interpretation": "bulk1_interpretation_db.json",
                    "population": "bulk1_population_db.json",
                    "epigenetics": "BULK1_epigenetics_hg19.csv",
                },
            },
            "CURATED1": {"shard": shard_name, "files": {"interpretation": "curated1_interpretation_db.json"}},
        },
        "files": {
            "bulk1_interpretation_db.json": shard_name,
            "bulk1_population_db.json": shard_name,
            "BULK1_epigenetics_hg19.csv": shard_name,
            "curated1_interpretation_db.json": shard_name,
        },
        "shards": {shard_name: {"path": f"bulk_gene_data_shards/{shard_name}"}},
    }
    analysis.GENE_DATA_INDEX_PATH.write_text(json.dumps(index), encoding="utf-8")
    analysis.clear_gene_data_bundle_cache()


def test_bulk_index_lookup_and_bulk_only_loaders(monkeypatch, tmp_path: Path) -> None:
    _install_bulk_test_data(monkeypatch, tmp_path)

    assert "bulk1_interpretation_db.json" in analysis.list_gene_data_bulk_members("_interpretation_db.json")

    knowledge_base = analysis.load_gene_interpretation_database("BULK1")
    population = analysis.load_gene_population_database("BULK1")
    manifest = analysis.load_gene_epigenetics_manifest("BULK1")

    assert knowledge_base is not None
    assert population is not None
    assert manifest is not None
    assert manifest.empty
    assert knowledge_base["gene_context"]["gene_region"]["start"] == 101
    assert population["database_name"] == "bulk population"


def test_curated_bundle_and_loose_files_override_bulk(monkeypatch, tmp_path: Path) -> None:
    _install_bulk_test_data(monkeypatch, tmp_path)

    curated = analysis.load_gene_interpretation_database("CURATED1")
    assert curated is not None
    assert curated["gene_context"]["gene_region"]["start"] == 11

    loose = _minimal_interpretation("BULK1", 303)
    (analysis.GENE_DATA_DIR / "bulk1_interpretation_db.json").write_text(json.dumps(loose), encoding="utf-8")
    overridden = analysis.load_gene_interpretation_database("BULK1")
    assert overridden is not None
    assert overridden["gene_context"]["gene_region"]["start"] == 303


def test_bulk_generator_is_deterministic_for_same_inputs(tmp_path: Path) -> None:
    hgnc_path = tmp_path / "hgnc.tsv"
    hgnc_path.write_text(
        "\t".join(["hgnc_id", "symbol", "name", "locus_group", "locus_type", "status", "location", "entrez_id", "ensembl_gene_id"])
        + "\n"
        + "HGNC:1\tALPHA1\talpha one\tprotein-coding gene\tgene with protein product\tApproved\t1p36\t1\tENSG000001\n"
        + "HGNC:2\tBETA1\tbeta one\tprotein-coding gene\tgene with protein product\tApproved\t2q11\t2\tENSG000002\n",
        encoding="utf-8",
    )
    gtf_path = tmp_path / "genes.gtf.gz"
    with gzip.open(gtf_path, "wt", encoding="utf-8") as handle:
        handle.write(
            'chr1\tTEST\tgene\t100\t200\t.\t+\t.\tgene_id "ENSG000001.1"; gene_type "protein_coding"; gene_name "ALPHA1"; gene_status "KNOWN";\n'
            'chr2\tTEST\tgene\t500\t800\t.\t-\t.\tgene_id "ENSG000002.1"; gene_type "protein_coding"; gene_name "BETA1"; gene_status "KNOWN";\n'
        )
    manifest_path = tmp_path / "manifest.csv"
    manifest_path.write_text(
        "IlmnID,CHR,MAPINFO,UCSC_RefGene_Name,UCSC_RefGene_Group,Relation_to_UCSC_CpG_Island,GencodeBasicV12_NAME,SNP_ID,SNP_DISTANCE\n"
        "cg00000001,1,150,ALPHA1,TSS200,Island,ALPHA1,rs1,0\n",
        encoding="utf-8",
    )

    first = build_bulk_gene_data(
        hgnc_path=hgnc_path,
        gencode_path=gtf_path,
        manifest_path=manifest_path,
        curated_bundle_path=tmp_path / "missing.zip",
        output_dir=tmp_path / "out1",
        shard_count=4,
    )
    second = build_bulk_gene_data(
        hgnc_path=hgnc_path,
        gencode_path=gtf_path,
        manifest_path=manifest_path,
        curated_bundle_path=tmp_path / "missing.zip",
        output_dir=tmp_path / "out2",
        shard_count=4,
    )

    assert first["files"] == second["files"]
    assert {name: data["sha256"] for name, data in first["shards"].items()} == {
        name: data["sha256"] for name, data in second["shards"].items()
    }


def test_actual_bulk_index_integrity_when_generated() -> None:
    index_path = analysis.GENE_DATA_INDEX_PATH
    assert index_path.exists()
    index = json.loads(index_path.read_text(encoding="utf-8"))
    assert index["bulk_gene_count"] > 10_000

    curated_genes: set[str] = set()
    with zipfile.ZipFile(analysis.GENE_DATA_BUNDLE_PATH) as curated_bundle:
        for member in curated_bundle.namelist():
            if member.endswith("_interpretation_db.json"):
                payload = json.loads(curated_bundle.read(member).decode("utf-8"))
                curated_genes.add(payload["gene_context"]["gene_name"].upper())

    for gene_name, entry in index["genes"].items():
        assert gene_name.upper() not in curated_genes
        assert set(entry["files"]) == {"interpretation", "population", "epigenetics"}
        for filename in entry["files"].values():
            assert index["files"][filename] == entry["shard"]

    required_columns = {"IlmnID", "CHR", "MAPINFO", "UCSC_RefGene_Name"}
    seen_files: set[str] = set()
    for shard_name, shard_info in index["shards"].items():
        shard_path = analysis.GENE_DATA_DIR / shard_info["path"]
        assert shard_path.exists(), shard_name
        with zipfile.ZipFile(shard_path) as bundle:
            members = {info.filename for info in bundle.infolist() if not info.is_dir()}
            seen_files.update(members)
            for member in members:
                if member.endswith(".json"):
                    json.loads(bundle.read(member).decode("utf-8"))
                elif member.endswith(".csv"):
                    header = pd.read_csv(bundle.open(member), nrows=0).columns
                    assert required_columns.issubset(set(header)), member
    assert seen_files == set(index["files"])
