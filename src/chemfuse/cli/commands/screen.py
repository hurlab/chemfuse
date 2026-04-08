"""Screen command for ChemFuse CLI.

Runs the full compound screening pipeline:
  read -> search -> enrich -> ADMET -> filter -> cluster -> export
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import click

from chemfuse.cli._utils import _validate_output_path
from chemfuse.models.collection import CompoundCollection

logger = logging.getLogger(__name__)


@click.command("screen")
@click.argument("input_file", type=click.Path(exists=True, readable=True))
@click.option(
    "--sources",
    default="pubchem",
    show_default=True,
    help="Comma-separated source databases (e.g. pubchem,chembl).",
)
@click.option(
    "--admet",
    "run_admet",
    is_flag=True,
    default=False,
    help="Run ADMET prediction on all compounds.",
)
@click.option(
    "--druglikeness",
    "druglikeness_rules",
    default=None,
    help="Comma-separated drug-likeness rules to apply (e.g. lipinski,veber).",
)
@click.option(
    "--cluster",
    "cluster_method",
    default=None,
    type=click.Choice(["butina", "kmeans"], case_sensitive=False),
    help="Clustering method to apply.",
)
@click.option(
    "--cluster-cutoff",
    default=0.4,
    show_default=True,
    type=float,
    help="Butina distance cutoff.",
)
@click.option(
    "--n-clusters",
    default=5,
    show_default=True,
    type=int,
    help="Number of clusters for KMeans.",
)
@click.option(
    "--output",
    "-o",
    default=None,
    help="Output file path (CSV, JSON, or XLSX).",
)
@click.option(
    "--limit",
    "-n",
    default=10,
    show_default=True,
    type=int,
    help="Maximum results per source per query.",
)
@click.pass_context
def screen_cmd(
    ctx: click.Context,
    input_file: str,
    sources: str,
    run_admet: bool,
    druglikeness_rules: str | None,
    cluster_method: str | None,
    cluster_cutoff: float,
    n_clusters: int,
    output: str | None,
    limit: int,
) -> None:
    """Run the full compound screening pipeline on INPUT_FILE.

    INPUT_FILE should be a CSV or text file with one SMILES or compound name
    per line (or a 'smiles' / 'name' column in CSV format).

    Examples:

        chemfuse screen compounds.csv --admet --cluster butina -o results.xlsx

        chemfuse screen smiles.txt --sources pubchem,chembl --druglikeness lipinski
    """
    # Rich progress import (graceful fallback)
    try:
        from rich.console import Console
        from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn

        console = Console()
        _rich = True
    except ImportError:
        console = None  # type: ignore[assignment]
        _rich = False

    def _echo(msg: str) -> None:
        if _rich and console is not None:
            console.print(msg)
        else:
            click.echo(msg, err=True)

    source_list = [s.strip() for s in sources.split(",") if s.strip()]

    # ------------------------------------------------------------------
    # Step 1: Read input file
    # ------------------------------------------------------------------
    _echo(f"[bold cyan]ChemFuse Screen[/bold cyan]: reading {input_file}" if _rich else f"Reading {input_file}...")
    queries = _read_input(input_file)
    if not queries:
        _echo("No compounds found in input file.")
        sys.exit(1)
    _echo(f"  Loaded {len(queries)} compound(s) from input.")

    # ------------------------------------------------------------------
    # Step 2: Batch search
    # ------------------------------------------------------------------
    all_compounds: list = []
    _echo(f"  Searching {len(queries)} queries across: {', '.join(source_list)}")
    if _rich and console is not None:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("{task.completed}/{task.total}"),
            console=console,
        ) as progress:
            task = progress.add_task("Searching sources...", total=len(queries))
            for query in queries:
                compounds = _search_single(query, source_list, limit)
                all_compounds.extend(compounds)
                progress.advance(task)
    else:
        for i, query in enumerate(queries, 1):
            compounds = _search_single(query, source_list, limit)
            all_compounds.extend(compounds)
            click.echo(f"  [{i}/{len(queries)}] {query}: {len(compounds)} result(s)", err=True)

    collection = CompoundCollection(compounds=all_compounds)
    _echo(f"  Total compounds after search: {len(collection)}")

    if not collection:
        _echo("No compounds retrieved.")
        sys.exit(0)

    # ------------------------------------------------------------------
    # Step 3: Enrich (compute descriptors + druglikeness)
    # ------------------------------------------------------------------
    _echo("  Computing descriptors...")
    _run_with_progress(
        collection.compute_all,
        description="Computing descriptors",
        use_rich=_rich,
        console=console,
    )

    # ------------------------------------------------------------------
    # Step 4: ADMET prediction
    # ------------------------------------------------------------------
    admet_count = 0
    if run_admet:
        _echo("  Predicting ADMET properties...")
        collection.predict_admet()
        admet_count = len(collection)
        _echo(f"  ADMET predicted for {admet_count} compound(s).")

    # ------------------------------------------------------------------
    # Step 5: Drug-likeness filter
    # ------------------------------------------------------------------
    if druglikeness_rules:
        rules = [r.strip().lower() for r in druglikeness_rules.split(",") if r.strip()]
        before = len(collection)
        for rule in rules:
            collection = collection.filter_druglike(rule=rule)
        after = len(collection)
        _echo(f"  Drug-likeness filter ({', '.join(rules)}): {before} -> {after} compound(s).")

    # ------------------------------------------------------------------
    # Step 6: Cluster
    # ------------------------------------------------------------------
    cluster_labels: list[int] = []
    n_clusters_found = 0
    if cluster_method and len(collection) >= 2:
        _echo(f"  Clustering with {cluster_method}...")
        try:
            cluster_labels = collection.cluster(
                method=cluster_method,
                cutoff=cluster_cutoff,
                n_clusters=n_clusters,
            )
            n_clusters_found = len(set(cluster_labels))
            _echo(f"  Found {n_clusters_found} cluster(s).")
        except Exception as exc:
            _echo(f"  Clustering failed (optional): {exc}")

    # ------------------------------------------------------------------
    # Step 7: Export
    # ------------------------------------------------------------------
    if output:
        _validate_output_path(Path(output))
        ext = Path(output).suffix.lower().lstrip(".")
        if ext == "json":
            collection.to_json(output)
        elif ext in ("xlsx", "xls"):
            collection.to_excel(output)
        else:
            collection.to_csv(output)
        _echo(f"  Results saved to: {output}")
    else:
        # Print brief summary table
        _print_summary(collection, use_rich=_rich, console=console)

    # ------------------------------------------------------------------
    # Summary statistics
    # ------------------------------------------------------------------
    _echo("")
    _echo("=== Screen Summary ===")
    _echo(f"  Input queries       : {len(queries)}")
    _echo(f"  Compounds found     : {len(collection)}")
    _echo(f"  ADMET predicted     : {admet_count}")
    _echo(f"  Clusters found      : {n_clusters_found or 'N/A'}")
    if output:
        _echo(f"  Output file         : {output}")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_input(path: str) -> list[str]:
    """Read compound queries from a file.

    Supports CSV files with 'smiles' or 'name' column, or plain text
    (one query per line).
    """
    import csv

    queries: list[str] = []
    p = Path(path)

    if p.suffix.lower() == ".csv":
        with p.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames:
                smiles_col = next(
                    (col for col in reader.fieldnames if col.lower() in ("smiles", "smi")), None
                )
                name_col = next(
                    (col for col in reader.fieldnames if col.lower() in ("name", "compound_name", "compound")), None
                )
                col = smiles_col or name_col
                if col:
                    for row in reader:
                        val = row.get(col, "").strip()
                        if val:
                            queries.append(val)
                    return queries
            # Fallback: treat first column as queries
            f.seek(0)
            for row in csv.reader(f):
                if row:
                    val = row[0].strip()
                    if val and not val.lower().startswith(("smiles", "name", "compound")):
                        queries.append(val)
    else:
        for line in p.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                queries.append(line)

    return queries


def _search_single(query: str, sources: list[str], limit: int) -> list:
    """Search a single query across sources. Returns list of Compound objects."""
    try:
        from chemfuse import search_async
        from chemfuse.cli._async import _run_async

        collection = _run_async(search_async(query, sources=sources, limit=limit))
        return list(collection.compounds) if collection else []
    except Exception as exc:
        logger.debug("Search failed for %r: %s", query, exc)
        return []


def _run_with_progress(
    func: object,
    description: str,
    use_rich: bool,
    console: object,
) -> None:
    """Run a callable with optional rich progress display."""
    if use_rich and console is not None:
        try:
            from rich.progress import Progress, SpinnerColumn, TextColumn

            with Progress(SpinnerColumn(), TextColumn(description), console=console) as p:
                p.add_task(description, total=None)
                func()  # type: ignore[operator]
            return
        except ImportError:
            pass
    func()  # type: ignore[operator]


def _print_summary(collection: object, use_rich: bool, console: object) -> None:
    """Print a brief summary of the collection."""
    if use_rich and console is not None:
        try:
            from rich.table import Table

            table = Table(title=f"Screen Results ({len(collection)} compounds)")  # type: ignore[arg-type]
            table.add_column("Name", style="cyan")
            table.add_column("SMILES", style="dim", max_width=40)
            table.add_column("MW", style="magenta")
            table.add_column("Sources", style="blue")

            for compound in collection:  # type: ignore[union-attr]
                mw = compound.properties.molecular_weight
                table.add_row(
                    compound.name or "",
                    (compound.smiles[:37] + "..." if len(compound.smiles) > 40 else compound.smiles),
                    f"{mw:.1f}" if mw else "",
                    ", ".join(compound.sources),
                )

            console.print(table)  # type: ignore[union-attr]
            return
        except ImportError:
            pass

    for i, compound in enumerate(collection, 1):  # type: ignore[call-overload]
        click.echo(f"{i:3}. {compound.name or 'Unknown':30} {compound.smiles[:50]}")
