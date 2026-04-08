"""CLI command to start the ChemFuse REST API server."""

from __future__ import annotations

import click


@click.command()
@click.option("--host", default="127.0.0.1", show_default=True, help="Host to bind to.")
@click.option("--port", default=8000, type=int, show_default=True, help="Port to bind to.")
@click.option("--reload", is_flag=True, help="Enable auto-reload for development.")
def serve(host: str, port: int, reload: bool) -> None:
    """Start the ChemFuse REST API server."""
    try:
        import uvicorn
    except ImportError as exc:
        click.echo(
            "Error: uvicorn is required. Install with: pip install chemfuse[api]",
            err=True,
        )
        raise SystemExit(1) from exc

    click.echo(f"Starting ChemFuse API at http://{host}:{port}")
    click.echo(f"API docs at http://{host}:{port}/docs")
    uvicorn.run(
        "chemfuse.api.server:create_app",
        host=host,
        port=port,
        reload=reload,
        factory=True,
    )
