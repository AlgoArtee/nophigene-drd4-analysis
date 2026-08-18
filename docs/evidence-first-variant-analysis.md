# Evidence-first variant analysis

NophiGene Version 2 treats a VCF row as an objective observation, not a phenotype statement.

## Required identity

Every variant retains its native assembly, chromosome, coordinate, REF, ALT, genotype, filtering fields, and raw source fields. Multi-allelic rows are split, normalized, left-aligned, and reference-validated before external adapters can consume them. GRCh37 and GRCh38 are different native identities.

Liftover records the original coordinate, target coordinate, chain, tool/version, ambiguity, status, and target REF validation. An ambiguous or REF-invalid mapping blocks downstream adapters.

## QC

Primary short-variant results require PASS plus GQ ≥20 and DP ≥10 when those fields are present. Reference calls and failed/flagged observations remain in the database, filters, Run Details, and full-fidelity exports.

SVs, CNVs, and VNTRs are reported as unsupported primary classes until a dedicated validated workflow exists.

## Evidence separation

- Objective Data contains the measured call.
- Statistics contains matched-cohort comparisons and uncertainty.
- Scientific Literature contains primary studies, reviews, preprints, trials, GWAS, and exploratory assertions.
- Medical Information accepts only authoritative reviewed clinical assertions with authority, jurisdiction, effective date/release, population, and applicability.
- Interactions contains typed source-backed gene edges.
- Predictions contains each computational model output independently.

No genotype, allele frequency, methylation mean, model score, or structure overlap is converted into a diagnosis. A prediction can appear in Medical only when a separate authoritative medical record supports the claim.

## Pharmacogenomics

Diplotypes are estimates in Predictions. Every defining locus must be QC-covered and exactly one solution must remain. Unphased, ambiguous, incomplete, or CNV-dependent calls remain unresolved, and phenotype is not assessed.
