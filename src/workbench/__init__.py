"""Evidence-first persistence, reporting, and adapter contracts."""

from .reporting import REPORT_SCHEMA_VERSION, build_canonical_report

__all__ = ["REPORT_SCHEMA_VERSION", "build_canonical_report"]
