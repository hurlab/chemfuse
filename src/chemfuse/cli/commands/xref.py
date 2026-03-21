"""Cross-reference command for ChemFuse CLI."""

from __future__ import annotations

import json
import sys

import click

from chemfuse.core.exceptions import ChemFuseError


@click.command("xref")
@click.option(
    "--cid",
    default=None,
    type=int,
    help="PubChem CID.",
)
@click.option(
    "--chembl",
    "chembl_id",
    default=None,
    help="ChEMBL ID (e.g., CHEMBL25).",
)
@click.option(
    "--smiles",
    default=None,
    help="SMILES string.",
)
@click.option(
    "--inchikey",
    default=None,
    help="Standard InChIKey (27 characters).",
)
@click.option(
    "--format",
    "fmt",
    default="table",
    show_default=True,
    type=click.Choice(["table", "json"], case_sensitive=False),
    help="Output format.",
)
@click.option(
    "--output",
    "-o",
    default=None,
    help="Output file path (JSON based on format).",
)
@click.pass_context
def xref_cmd(
    ctx: click.Context,
    cid: int | None,
    chembl_id: str | None,
    smiles: str | None,
    inchikey: str | None,
    fmt: str,
    output: str | None,
) -> None:
    """Cross-reference a compound identifier across chemical databases.

    At least one of --cid, --chembl, --smiles, or --inchikey must be provided.
    Uses UniChem to map identifiers between PubChem, ChEMBL, DrugBank, KEGG,
    and ChEBI.

    Examples:

        chemfuse xref --cid 2244

        chemfuse xref --chembl CHEMBL25

        chemfuse xref --inchikey BSYNRYMUTXBXSQ-UHFFFAOYSA-N

        chemfuse xref --cid 2244 --format json
    """
    if all(v is None for v in [cid, chembl_id, smiles, inchikey]):
        click.echo(
            "Error: At least one of --cid, --chembl, --smiles, or --inchikey is required.",
            err=True,
        )
        sys.exit(1)

    try:
        import asyncio

        from chemfuse import map_identifiers_async

        mapping = asyncio.run(
            map_identifiers_async(
                cid=cid,
                chembl_id=chembl_id,
                inchikey=inchikey,
                smiles=smiles,
            )
        )
    except ChemFuseError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"Unexpected error: {e}", err=True)
        sys.exit(1)

    if not mapping:
        click.echo("No cross-references found.", err=True)
        sys.exit(0)

    if output:
        with open(output, "w") as f:
            json.dump(mapping, f, indent=2)
        click.echo(f"Results saved to: {output}", err=True)
        return

    if fmt == "json":
        click.echo(json.dumps(mapping, indent=2))
    else:
        # Table format
        _print_xref_table(mapping)


def _print_xref_table(mapping: dict[str, str]) -> None:
    """Print cross-reference mapping as a formatted table."""
    try:
        from rich.console import Console
        from rich.table import Table

        console = Console()
        table = Table(title="Cross-Reference Mapping")
        table.add_column("Database", style="cyan", no_wrap=True)
        table.add_column("Identifier", style="green")

        for db_name, identifier in sorted(mapping.items()):
            table.add_row(db_name, identifier)

        console.print(table)
    except ImportError:
        # Fallback to simple text output
        click.echo("Database       Identifier")
        click.echo("-" * 40)
        for db_name, identifier in sorted(mapping.items()):
            click.echo(f"{db_name:<15} {identifier}")
