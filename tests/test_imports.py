"""Import smoke tests for the DRD4 analysis project."""

from src import analysis, app, gene_region_extraction, webapp


def test_analysis_module_imports() -> None:
    """Ensure the core analysis module remains importable."""
    assert analysis.DEFAULT_REGION


def test_web_modules_import() -> None:
    """Ensure the UI entrypoints remain importable."""
    assert app.build_parser() is not None
    assert webapp.app is not None


def test_gene_region_helpers_import() -> None:
    """Ensure the region helper module remains importable."""
    assert callable(gene_region_extraction.get_widest_region)
