# NophiGene Version 1 archive

This directory contains retired Version 1 material retained for historical reference and migration work. Nothing here is imported, copied into the application image, or used by the Version 2 launchers.

Tracked historical files are preserved byte-for-byte, so paths written inside archived launchers and documentation describe their original locations. Active Version 2 references were updated instead of rewriting the archive.

Archived material includes:

- the original monolithic result template and genotype-aware design note;
- the retired predictive-synthesis generator;
- the former local/Docker launcher variants and development-container configuration;
- DRD4 exploratory notebook/helper material;
- historical research notes, changelog, and optional research dependencies;
- loose per-gene JSON/CSV stores superseded by the active normalized bundles and Version 2 database.

The active Version 2 application remains at the repository root. Use `Start NophiGene UI.cmd`, `Stop NophiGene UI.cmd`, `scripts/start-v2.ps1`, and `scripts/stop-v2.ps1`.

## Deliberately retained outside this archive

Some code originated before Version 2 but is not a relic:

- `src/analysis.py`, preprocessing helpers, and the variant-knowledge adapters are still used by the Version 2 workflow while normalized services replace them incrementally.
- `src/templates/functional_map.html` remains the source-backed functional-family explorer required by Data Explorer; it is not presented as an interaction network.
- `src/api/routes.py` is the intentionally read-only `/api/v1` compatibility surface for legacy jobs and artifacts.
- the active gene-data bundle and sharded index remain required for reference lookup and migration verification.
- root `resources.txt` remains the active input catalog parsed by the Version 2 source registry.

Do not run archived launchers against sensitive samples. They are preserved to explain historical behavior, not as a supported runtime.
