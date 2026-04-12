"""Container entrypoint that launches either the CLI workflow or the web UI."""

from __future__ import annotations

import argparse

try:
    from . import analysis
    from .webapp import run_server
except ImportError:
    import analysis
    from webapp import run_server


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level launcher parser."""
    parser = argparse.ArgumentParser(
        description="Launch the NophiGene DRD4 analysis app in CLI or web mode."
    )
    subparsers = parser.add_subparsers(dest="mode")

    web_parser = subparsers.add_parser("web", help="Run the browser-based interface.")
    web_parser.add_argument("--host", default="0.0.0.0", help="Host interface to bind.")
    web_parser.add_argument("--port", default=8000, type=int, help="Port to serve the UI on.")
    web_parser.add_argument("--debug", action="store_true", help="Enable Flask debug mode.")

    cli_parser = subparsers.add_parser("cli", help="Run the original command-line pipeline.")
    cli_parser.add_argument("analysis_args", nargs=argparse.REMAINDER, help="Arguments forwarded to src/analysis.py")

    return parser


def main(argv: list[str] | None = None) -> int:
    """Dispatch to web mode or CLI mode based on the selected subcommand."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.mode in (None, "web"):
        host = getattr(args, "host", "0.0.0.0")
        port = getattr(args, "port", 8000)
        debug = getattr(args, "debug", False)
        run_server(host=host, port=port, debug=debug)
        return 0

    if args.mode == "cli":
        analysis_args = args.analysis_args
        if analysis_args and analysis_args[0] == "--":
            analysis_args = analysis_args[1:]
        return analysis.main(analysis_args)

    parser.error(f"Unsupported mode: {args.mode}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
