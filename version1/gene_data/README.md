# Retired loose gene stores

The local `loose/` subdirectory contains the Version 1 per-gene JSON/CSV files that formerly overrode the bundled store. It is intentionally Git-ignored because it contains thousands of generated artifacts. The small tracked sidecar files are under `metadata/`.

Relocated inventory:

- 637 synthetic synthesis JSON files (archive-only; never loaded by Version 2);
- 635 interpretation JSON files;
- 635 population JSON files;
- 636 epigenetics CSV subsets;
- 20 sidecar metadata files.

Every relocated interpretation, population, and epigenetics file had a corresponding member in `src/gene_data/gene_data_bundle.zip` before relocation. Version 2 therefore continues to resolve those records from the active bundle, while synthetic synthesis records have no active loader.
