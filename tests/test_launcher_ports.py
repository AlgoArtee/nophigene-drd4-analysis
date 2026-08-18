"""Regression checks for consistent UI port configuration."""

from pathlib import Path

from src import app, webapp

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WEB_PORT = 8766


def test_python_web_defaults_use_shared_preferred_port() -> None:
    args = app.build_parser().parse_args(["web"])

    assert app.DEFAULT_WEB_PORT == DEFAULT_WEB_PORT
    assert args.port == DEFAULT_WEB_PORT
    assert webapp.run_server.__defaults__ == ("127.0.0.1", DEFAULT_WEB_PORT, False)


def test_launchers_are_verbose_version_two_compose_entry_points() -> None:
    start_script = (PROJECT_ROOT / "scripts" / "start-v2.ps1").read_text(encoding="utf-8")
    stop_script = (PROJECT_ROOT / "scripts" / "stop-v2.ps1").read_text(encoding="utf-8")

    assert "[int]$Port = 8766" in start_script
    assert "/api/v2/health" in start_script
    assert "Write-Stage" in start_script
    assert "Show-ComposeDiagnostics" in start_script
    assert '"config", "--quiet"' in start_script
    assert '"Secret values" "redacted"' in start_script
    assert "RandomNumberGenerator]::Create()" in start_script
    assert "RandomNumberGenerator]::Fill" not in start_script
    assert "docker compose" in stop_script
    assert "Runtime secret" in stop_script or "runtime secret" in stop_script
    assert "Results retained" in stop_script

    for launcher_name in ("Start NophiGene UI.cmd", "Stop NophiGene UI.cmd"):
        launcher = (PROJECT_ROOT / launcher_name).read_text(encoding="utf-8")
        assert "Version 2" in launcher
        assert "Exit code" in launcher
        assert "TARGET_SCRIPT" in launcher

    assert not (PROJECT_ROOT / "Start NophiGene UI (Docker).cmd").exists()
    assert not (PROJECT_ROOT / "Stop NophiGene UI (Docker).cmd").exists()
    assert (PROJECT_ROOT / "version1" / "launchers" / "Start NophiGene UI (Docker).cmd").exists()
    assert (PROJECT_ROOT / "version1" / "launchers" / "Stop NophiGene UI (Docker).cmd").exists()
    assert (PROJECT_ROOT / "version1" / "launchers" / "scripts" / "start_nophigene_ui_local.ps1").exists()
    assert (PROJECT_ROOT / "version1" / "launchers" / "scripts" / "stop_nophigene_ui_local.ps1").exists()


def test_container_and_ui_copy_do_not_reference_old_port() -> None:
    dockerfile = (PROJECT_ROOT / "Dockerfile").read_text(encoding="utf-8")
    template = (PROJECT_ROOT / "src" / "templates" / "v2" / "index.html").read_text(
        encoding="utf-8"
    )
    compose = (PROJECT_ROOT / "docker-compose.yml").read_text(encoding="utf-8")

    assert "EXPOSE 8766" in dockerfile
    assert '"--port", "8766"' in dockerfile
    assert '127.0.0.1:${NOPHIGENE_PORT:-8766}:8766' in compose
    assert "127.0.0.1:8000" not in template
