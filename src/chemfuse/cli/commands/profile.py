"""Profile command for ChemFuse CLI - detailed compound information."""

from __future__ import annotations

import sys

import click

from chemfuse.core.exceptions import NotFoundError, SourceError


@click.command("profile")
@click.argument("identifier")
@click.option(
    "--source",
    "-s",
    default="pubchem",
    show_default=True,
    help="Source database.",
)
@click.option(
    "--similar",
    is_flag=True,
    default=False,
    help="Also find similar compounds.",
)
@click.option(
    "--threshold",
    default=90,
    show_default=True,
    type=int,
    help="Similarity threshold (0-100) for --similar.",
)
@click.option(
    "--patents",
    is_flag=True,
    default=False,
    help="Include patent data from SureChEMBL.",
)
@click.option(
    "--targets",
    is_flag=True,
    default=False,
    help="Include disease-target associations from Open Targets.",
)
@click.option(
    "--binding",
    is_flag=True,
    default=False,
    help="Include binding data from BindingDB.",
)
@click.option(
    "--all",
    "all_sources",
    is_flag=True,
    default=False,
    help="Enable all enrichment sources (equivalent to --patents --targets --binding).",
)
@click.pass_context
def profile_cmd(
    ctx: click.Context,
    identifier: str,
    source: str,
    similar: bool,
    threshold: int,
    patents: bool,
    targets: bool,
    binding: bool,
    all_sources: bool,
) -> None:
    """Show detailed profile for a compound.

    IDENTIFIER is a database-specific compound identifier (e.g., PubChem CID).

    Examples:

        chemfuse profile 2244

        chemfuse profile 2244 --similar

        chemfuse profile 2244 --similar --threshold 85

        chemfuse profile aspirin --patents --targets

        chemfuse profile aspirin --all
    """
    # --all enables all enrichment sources
    if all_sources:
        patents = True
        targets = True
        binding = True

    try:
        from chemfuse import get_async
        from chemfuse.cli._async import _run_async

        compound = _run_async(get_async(identifier, source=source))
    except NotFoundError as e:
        click.echo(f"Not found: {e}", err=True)
        sys.exit(1)
    except SourceError as e:
        click.echo(f"Source error: {e}", err=True)
        sys.exit(1)

    if compound is None:
        click.echo(f"Compound '{identifier}' not found in {source}.", err=True)
        sys.exit(1)

    # Perform enrichment if requested
    if patents or targets or binding:
        try:
            from chemfuse.cli._async import _run_async
            _run_async(compound.enrich(patents=patents, targets=targets, binding=binding))
        except Exception as e:
            click.echo(f"Enrichment warning: {e}", err=True)

    _print_profile(compound)

    if similar and compound.smiles:
        click.echo("\nFinding similar compounds...")
        try:
            from chemfuse import find_similar_async
            from chemfuse.cli._async import _run_async

            similar_collection = _run_async(
                find_similar_async(
                    compound.smiles,
                    threshold=threshold,
                    max_results=10,
                    sources=[source],
                )
            )
            if similar_collection:
                click.echo(f"\nSimilar compounds (threshold={threshold}%):")
                for i, sim_compound in enumerate(similar_collection, 1):
                    name = sim_compound.name or "Unknown"
                    cid = f"CID:{sim_compound.cid}" if sim_compound.cid else ""
                    click.echo(f"  {i:2}. {name:30} {cid}")
            else:
                click.echo("No similar compounds found.")
        except Exception as e:
            click.echo(f"Similarity search failed: {e}", err=True)


def _print_profile(compound: object) -> None:
    """Print compound profile."""
    try:
        from rich import box
        from rich.console import Console
        from rich.panel import Panel
        from rich.table import Table

        console = Console()

        name = getattr(compound, "name", None) or "Unknown"
        panel = Panel(
            f"[bold cyan]{name}[/bold cyan]",
            title="Compound Profile",
            border_style="blue",
        )
        console.print(panel)

        # Identifiers table
        id_table = Table(title="Identifiers", box=box.SIMPLE)
        id_table.add_column("Field", style="green")
        id_table.add_column("Value", style="white")

        fields = [
            ("CID", getattr(compound, "cid", None)),
            ("ChEMBL ID", getattr(compound, "chembl_id", None)),
            ("Formula", getattr(compound, "formula", None)),
            ("InChIKey", getattr(compound, "inchikey", None)),
            ("SMILES", getattr(compound, "smiles", None)),
        ]
        for field_name, value in fields:
            if value:
                id_table.add_row(field_name, str(value))
        console.print(id_table)

        # Properties table
        props = getattr(compound, "properties", None)
        if props:
            prop_table = Table(title="Properties", box=box.SIMPLE)
            prop_table.add_column("Property", style="green")
            prop_table.add_column("Value", style="white")

            prop_fields = [
                ("Molecular Weight", getattr(props, "molecular_weight", None), ".2f"),
                ("Exact Mass", getattr(props, "exact_mass", None), ".4f"),
                ("XLogP", getattr(props, "xlogp", None), ".2f"),
                ("TPSA", getattr(props, "tpsa", None), ".1f"),
                ("HBD Count", getattr(props, "hbd_count", None), "d"),
                ("HBA Count", getattr(props, "hba_count", None), "d"),
                ("Rotatable Bonds", getattr(props, "rotatable_bonds", None), "d"),
                ("Heavy Atom Count", getattr(props, "heavy_atom_count", None), "d"),
                ("Complexity", getattr(props, "complexity", None), ".1f"),
            ]
            for prop_name, value, fmt in prop_fields:
                if value is not None:
                    formatted = f"{value:{fmt}}"
                    prop_table.add_row(prop_name, formatted)
            console.print(prop_table)

        # Patents section
        patents_data = getattr(compound, "patents", [])
        if patents_data:
            pat_table = Table(title="Patents", box=box.SIMPLE)
            pat_table.add_column("Patent ID", style="green")
            pat_table.add_column("Title", style="white")
            pat_table.add_column("Filing Date", style="cyan")
            for patent in patents_data:
                pat_table.add_row(
                    getattr(patent, "patent_id", "") or "",
                    (getattr(patent, "title", "") or "")[:60],
                    getattr(patent, "filing_date", "") or "",
                )
            console.print(pat_table)

        # Disease Associations section
        target_assocs = getattr(compound, "target_associations", [])
        if target_assocs:
            ta_table = Table(title="Disease Associations", box=box.SIMPLE)
            ta_table.add_column("Target", style="green")
            ta_table.add_column("Disease", style="white")
            ta_table.add_column("Score", style="cyan")
            for assoc in target_assocs:
                score_str = ""
                score = getattr(assoc, "association_score", None)
                if score is not None:
                    score_str = f"{score:.3f}"
                ta_table.add_row(
                    getattr(assoc, "target_name", "") or "",
                    getattr(assoc, "disease_name", "") or "",
                    score_str,
                )
            console.print(ta_table)

        # Sources
        sources = getattr(compound, "sources", [])
        if sources:
            console.print(f"[dim]Sources: {', '.join(sources)}[/dim]")

    except ImportError:
        # Fallback plain text
        _print_profile_plain(compound)


def _print_profile_plain(compound: object) -> None:
    """Print compound profile as plain text."""
    click.echo("=" * 60)
    click.echo(f"Name: {getattr(compound, 'name', 'Unknown')}")
    click.echo("=" * 60)

    if cid := getattr(compound, "cid", None):
        click.echo(f"CID: {cid}")
    if formula := getattr(compound, "formula", None):
        click.echo(f"Formula: {formula}")
    if inchikey := getattr(compound, "inchikey", None):
        click.echo(f"InChIKey: {inchikey}")
    if smiles := getattr(compound, "smiles", None):
        click.echo(f"SMILES: {smiles}")

    props = getattr(compound, "properties", None)
    if props:
        click.echo("\nProperties:")
        if (mw := getattr(props, "molecular_weight", None)) is not None:
            click.echo(f"  Molecular Weight: {mw:.2f}")
        if (xlogp := getattr(props, "xlogp", None)) is not None:
            click.echo(f"  XLogP: {xlogp:.2f}")
        if (tpsa := getattr(props, "tpsa", None)) is not None:
            click.echo(f"  TPSA: {tpsa:.1f}")
        if (hbd := getattr(props, "hbd_count", None)) is not None:
            click.echo(f"  HBD Count: {hbd}")
        if (hba := getattr(props, "hba_count", None)) is not None:
            click.echo(f"  HBA Count: {hba}")

    sources = getattr(compound, "sources", [])
    if sources:
        click.echo(f"\nSources: {', '.join(sources)}")

    patents_data = getattr(compound, "patents", [])
    if patents_data:
        click.echo("\nPatents:")
        for patent in patents_data:
            pid = getattr(patent, "patent_id", "") or ""
            title = (getattr(patent, "title", "") or "")[:60]
            fdate = getattr(patent, "filing_date", "") or ""
            click.echo(f"  {pid}: {title} ({fdate})")

    target_assocs = getattr(compound, "target_associations", [])
    if target_assocs:
        click.echo("\nDisease Associations:")
        for assoc in target_assocs:
            tname = getattr(assoc, "target_name", "") or ""
            dname = getattr(assoc, "disease_name", "") or ""
            score = getattr(assoc, "association_score", None)
            score_str = f"{score:.3f}" if score is not None else "N/A"
            click.echo(f"  {tname} -> {dname}: {score_str}")
