# DANDELION integration

NophiGene runs [DANDELION 0.1.0](https://CRAN.R-project.org/package=DANDELION)
as an offline, cohort-level statistical workflow. DANDELION integrates
trans-regulatory association p-values with named gene-level disease-association
p-values. Its output is hypothesis-generating and cannot by itself establish
causality, pathogenicity, diagnosis, or treatment relevance.

## Runtime boundary

- R 4.6.1 and DANDELION 0.1.0 run in `dandelion-worker`.
- The container has no network, Docker socket, or Linux capabilities.
- The app submits HMAC-SHA256 authenticated JSON envelopes. The signing key is
  retrieved from Windows Credential Manager and mounted as a runtime secret.
- The worker re-hashes every input before running R. It accepts only relative
  paths below the read-only `data/dandelion` import root.
- One CPU job runs at a time. Cancellation is checked between chunks and while
  the R process is active; jobs have a configurable timeout.
- Delimited matrices are read in exposure chunks. RDS matrices must be loaded as
  one object, so their preflight estimate carries an explicit memory warning.
- The worker does not install `qvalue`; DANDELION therefore uses its documented
  Benjamini–Hochberg fallback. The runner recomputes BH q-values and aborts if
  package significance indicators do not match.

## Dataset manifest

Paths are relative to `data/dandelion`. The three required roles for gene
exposures are `trans_matrix`, `gene_association`, and `gene_annotation`. SNP
exposures also require `snp_reference`; `snp_gene_map` and `exposure_list` are
optional.

```json
{
  "name": "Example cohort",
  "description": "Trans statistics plus rare-variant burden results",
  "phenotype": "EFO:0000270",
  "exposure_type": "Gene",
  "assembly": "GRCh38",
  "gene_namespace": "HGNC symbol",
  "context": {
    "tissue": "whole blood",
    "ancestry": "declared cohort ancestry",
    "assay": "trans-eQTL summary statistics"
  },
  "files": {
    "trans_matrix": {
      "path": "example/trans.tsv",
      "format": "tsv",
      "mapping": {"row_id": "gene"}
    },
    "gene_association": {
      "path": "example/gene-burden.tsv",
      "mapping": {"gene": "gene", "p_value": "p_value"}
    },
    "gene_annotation": {
      "path": "example/gene-annotation.tsv",
      "mapping": {
        "gene_name": "gene_name",
        "type": "type",
        "Chromosome": "Chromosome",
        "start": "start",
        "end": "end"
      }
    }
  }
}
```

The trans matrix contains candidate disease-proximal genes in rows and Gene or
SNP exposures in columns. All p-values must be finite values in `[0, 1]` or
missing. Gene identifiers must match among the trans-matrix row names, the
gene-association file, and the annotation file. The annotation `type` must use
the values expected by DANDELION, including `protein_coding` or `lincRNA`, and
chromosomes are normalized to `chr1`-style values.

For an SNP analysis, `snp_reference` maps `SNP`, `SNPPos`, and `SNPChr`. An
optional `snp_gene_map` maps `SNP` to `GeneSymbol` for gene-first interaction
display. Without that optional mapping, NophiGene preserves the variant as the
source node and labels its node type rather than inventing a gene mapping.

## API

Register a manifest:

```http
POST /api/v2/statistical-datasets
Content-Type: application/json
```

Queue an analysis:

```json
{
  "method": "dandelion",
  "dataset_id": "DATASET_ID",
  "parameters": {
    "target_fdr": 0.1,
    "cis_window_bp": 5000000,
    "gene_association_threshold": 0.00001,
    "chunk_size": 250
  }
}
```

Submit that object to `POST /api/v2/statistical-analyses`. Poll
`GET /api/v2/statistical-analyses/{id}` and retrieve canonical report schema
3.0 from `GET /api/v2/statistical-analyses/{id}/result`. Cancellation uses
`POST /api/v2/statistical-analyses/{id}/cancel`.

`POST /api/v2/statistical-datasets/{id}/managed-copy` creates a new immutable
catalog entry backed by explicit AES-256 encrypted copies. The source catalog
entry and original files are left untouched.

## Result semantics

- Statistics includes the original trans p-value, gene-level disease p-value,
  DANDELION p-value, BH q-value, declared family, context, and limitations.
- Interactions includes directed trans-regulatory disease-prioritisation edges.
  The source-native score is the BH q-value and is labeled “lower is stronger”;
  it is never converted into a probability.
- Up to 20 non-significant ranked pairs per exposure are normalized for
  interactive review. All significant pairs and complete per-chunk package RDS
  objects are retained as artifacts.
- No asthma-specific example or claim is bundled. Asthma results appear only
  when the user registers a verified asthma dataset and declares its phenotype.
