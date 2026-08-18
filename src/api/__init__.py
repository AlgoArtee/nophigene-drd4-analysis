"""Versioned local REST API for NophiGene."""

from __future__ import annotations

from flask import Flask

from .jobs import get_default_job_manager
from .profiles import get_default_profile_store
from .routes import api_v1
from .v2_routes import api_v2
try:
    from ..workbench.database import get_default_engine
    from ..workbench.model_jobs import get_default_model_job_manager
    from ..workbench.statistical_jobs import get_default_statistical_job_manager
except ImportError:
    from workbench.database import get_default_engine
    from workbench.model_jobs import get_default_model_job_manager
    from workbench.statistical_jobs import get_default_statistical_job_manager


def register_api(app: Flask) -> None:
    """Register the API blueprint and start its local background worker."""
    if "api_v1" in app.blueprints and "api_v2" in app.blueprints:
        return
    app.config.setdefault("NOPHIGENE_PROFILE_STORE", get_default_profile_store())
    app.config.setdefault("NOPHIGENE_JOB_MANAGER", get_default_job_manager())
    # Open the encrypted database during startup while the mounted Compose
    # secret is guaranteed to exist.  This also fails closed before serving a
    # request if SQLCipher or the key is unavailable.
    app.config.setdefault("NOPHIGENE_DATABASE_ENGINE", get_default_engine())
    app.config.setdefault("NOPHIGENE_MODEL_JOB_MANAGER", get_default_model_job_manager())
    app.config.setdefault("NOPHIGENE_STATISTICAL_JOB_MANAGER", get_default_statistical_job_manager())
    app.register_blueprint(api_v1)
    app.register_blueprint(api_v2)
    app.config["NOPHIGENE_JOB_MANAGER"].start()
