"""CLI command to launch the ChemFuse Streamlit web UI."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import click

from chemfuse.core.exceptions import OptionalDependencyError


@click.command("web")
@click.option(
    "--port",
    "-p",
    default=int(os.environ.get("CHEMFUSE_WEB_PORT", "8501")),
    show_default=True,
    type=int,
    help="Port for the Streamlit server.",
)
@click.option(
    "--host",
    default="localhost",
    show_default=True,
    type=str,
    help="Host/address for the Streamlit server.",
)
@click.option(
    "--no-browser",
    is_flag=True,
    default=False,
    help="Do not open the browser automatically.",
)
def web_cmd(port: int, host: str, no_browser: bool) -> None:
    """Launch the ChemFuse web dashboard (Streamlit).

    Starts a Streamlit server and opens the ChemFuse web UI in the
    default browser.

    Requires the [web] optional dependency:
        pip install chemfuse[web]
    """
    # Check if streamlit is installed
    try:
        import streamlit  # noqa: F401
    except ImportError as exc:
        raise OptionalDependencyError("streamlit", extra="web") from exc

    # Locate app.py
    app_file = Path(__file__).parent.parent.parent / "web" / "app.py"

    # Build the streamlit run command
    cmd = [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        str(app_file),
        "--server.port",
        str(port),
        "--server.address",
        host,
    ]

    if no_browser:
        cmd += ["--server.headless", "true"]

    click.echo(f"Starting ChemFuse web UI on http://{host}:{port}")
    click.echo("Press Ctrl+C to stop the server.")

    try:
        subprocess.run(cmd, check=True)
    except KeyboardInterrupt:
        click.echo("\nChemFuse web server stopped.")
    except subprocess.CalledProcessError as exc:
        click.echo(f"Streamlit server exited with error: {exc}", err=True)
        sys.exit(1)
