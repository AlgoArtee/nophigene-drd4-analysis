# Changelog

## Version 2 — Unreleased

- Replaced the monolithic result with canonical report schema 3.0 and six evidence-separated result sections.
- Added normalized SQLAlchemy/Alembic persistence, SQLCipher container enforcement, audit chaining, backups, exports, and migration tooling.
- Added `/api/v2`; retained `/api/v1` as a read-only legacy compatibility surface.
- Removed Predictive Theses, synthetic interpretation cases, arbitrary methylation phenotype buckets, and hard-coded HERC2 phenotype rules.
- Added typed literature, medical-evidence gates, statistics safeguards, interaction graphs, PGx resolution, and independent model contracts.
- Consolidated supported execution on the secured Version 2 Docker Compose launcher and Python 3.12 `.venv-v2` development environment.
- Moved retired Version 1 material into `version1/`.

The historical Version 1 changelog is preserved at `version1/CHANGELOG.md`.
