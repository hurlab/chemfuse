"""Search command for ChemFuse CLI."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click

from chemfuse.core.exceptions import NotFoundError, SourceError


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


@click.command("search")
@click.argument("query")
@click.option(
    "--type",
    "query_type",
    default="name",
    show_default=True,
    type=click.Choice(["name", "smiles", "cid", "formula", "inchi"], case_sensitive=False),
    help="Query type.",
)
@click.option(
    "--source",
    "-s",
    "sources",
    multiple=True,
    default=["pubchem"],
    show_default=True,
    help="Source database(s) to search. Can be specified multiple times.",
)
@click.option(
    "--sources",
    "sources_csv",
    default=None,
    help="Comma-separated source names (e.g., pubchem,chembl). Overrides --source.",
)
@click.option(
    "--limit",
    "-n",
    default=10,
    show_default=True,
    type=int,
    help="Maximum results per source.",
)
@click.option(
    "--output",
    "-o",
    default=None,
    help="Output file path (CSV, JSON, or XLSX based on extension).",
)
@click.option(
    "--format",
    "fmt",
    default="table",
    show_default=True,
    type=click.Choice(["table", "json", "csv"], case_sensitive=False),
    help="Output format for terminal display.",
)
@click.pass_context
def search_cmd(
    ctx: click.Context,
    query: str,
    query_type: str,
    sources: tuple[str, ...],
    sources_csv: str | None,
    limit: int,
    output: str | None,
    fmt: str,
) -> None:
    """Search chemical databases for a compound.

    QUERY can be a compound name, SMILES, CID, molecular formula, or InChI.

    Examples:

        chemfuse search aspirin

        chemfuse search "CC(=O)Oc1ccccc1C(=O)O" --type smiles

        chemfuse search 2244 --type cid

        chemfuse search C9H8O4 --type formula
    """
    # Resolve sources: --sources CSV overrides --source flag
    effective_sources: list[str]
    if sources_csv is not None:
        effective_sources = [s.strip() for s in sources_csv.split(",") if s.strip()]
    else:
        effective_sources = list(sources)

    try:
        from chemfuse import search_async
        from chemfuse.cli._async import _run_async

        collection = _run_async(
            search_async(
                query,
                sources=effective_sources,
                query_type=query_type,
                limit=limit,
            )
        )
    except NotFoundError as e:
        click.echo(f"Not found: {e}", err=True)
        sys.exit(1)
    except SourceError as e:
        click.echo(f"Source error: {e}", err=True)
        sys.exit(1)

    if not collection:
        click.echo("No results found.", err=True)
        sys.exit(0)

    click.echo(f"Found {len(collection)} compound(s).", err=True)

    # Save to file if requested
    if output:
        _validate_output_path(Path(output))
        ext = output.rsplit(".", 1)[-1].lower() if "." in output else "csv"
        if ext == "json":
            collection.to_json(output)
        elif ext in ("xlsx", "xls"):
            collection.to_excel(output)
        else:
            collection.to_csv(output)
        click.echo(f"Results saved to: {output}", err=True)
        return

    # Terminal output
    if fmt == "json":
        records = []
        for compound in collection:
            records.append(compound.to_dict())
        click.echo(json.dumps(records, indent=2, default=str))
    elif fmt == "csv":
        df = collection.to_dataframe()
        click.echo(df.to_csv(index=False))
    else:
        # Table format using simple formatting
        _print_table(collection)


def _print_table(collection: object) -> None:
    """Print compounds as a simple text table."""
    try:
        from rich.console import Console
        from rich.table import Table

        console = Console()
        table = Table(title=f"ChemFuse Results ({len(collection)} compounds)")  # type: ignore[arg-type]
        table.add_column("Name", style="cyan", no_wrap=True)
        table.add_column("CID", style="green")
        table.add_column("Formula", style="yellow")
        table.add_column("MW", style="magenta")
        table.add_column("SMILES", style="dim", max_width=40)
        table.add_column("Sources", style="blue")

        for compound in collection:  # type: ignore[union-attr]
            table.add_row(
                compound.name or "",
                str(compound.cid) if compound.cid else "",
                compound.formula or "",
                f"{compound.properties.molecular_weight:.2f}" if compound.properties.molecular_weight else "",
                compound.smiles[:40] + "..." if len(compound.smiles) > 40 else compound.smiles,
                ", ".join(compound.sources),
            )

        console.print(table)
    except ImportError:
        # Fallback to simple text output
        for i, compound in enumerate(collection, 1):  # type: ignore[call-overload]
            name = compound.name or "Unknown"
            cid = f"CID:{compound.cid}" if compound.cid else ""
            mw = f"MW:{compound.properties.molecular_weight:.1f}" if compound.properties.molecular_weight else ""
            click.echo(f"{i:3}. {name:30} {cid:15} {mw:12} {compound.smiles[:50]}")
