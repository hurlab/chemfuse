"""Shared utilities for ChemFuse CLI commands."""

from __future__ import annotations

from pathlib import Path

import click


def _validate_output_path(path: Path) -> Path:
    """Validate output path is safe.

    Args:
        path: Output file path.

    Returns:
        Resolved path.

    Raises:
        click.BadParameter: If the path contains traversal components.
    """
    resolved = path.resolve()
    if ".." in path.parts:
        raise click.BadParameter(f"Path traversal not allowed: {path}")
    return resolved
