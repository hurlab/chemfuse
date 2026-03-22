"""Batch processing engine for compound searches from files."""

from __future__ import annotations

import asyncio
import csv
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from chemfuse.models.collection import CompoundCollection
    from chemfuse.models.compound import Compound


MAX_BATCH_SIZE = 10_000


def _detect_query_type(query: str) -> str:
    """Auto-detect query type from the query string.

    Args:
        query: The query string.

    Returns:
        One of 'cid', 'smiles', 'formula', or 'name'.
    """
    query = query.strip()

    # Pure digits -> CID
    if query.isdigit():
        return "cid"

    # InChI prefix
    if query.startswith("InChI="):
        return "inchi"

    # Contains typical SMILES characters -> SMILES
    smiles_indicators = {"(", ")", "=", "#", "[", "]", "/", "\\", "@", "+"}
    has_smiles_chars = any(c in smiles_indicators for c in query)
    has_organic_subset = any(
        pat in query for pat in ["C(", "c1", "C=", "O)", "N(", "CC", "c(", "C#"]
    )

    if has_smiles_chars or has_organic_subset:
        return "smiles"

    # Molecular formula pattern: starts with C followed by digits and element symbols
    import re
    formula_pattern = r"^[A-Z][a-z]?(\d+)?([A-Z][a-z]?(\d+)?)+$"
    if re.match(formula_pattern, query) and query[0].isupper():
        return "formula"

    # Default to name search
    return "name"


def _read_queries_from_file(input_path: str | Path) -> list[str]:
    """Read compound queries from a file.

    Supports TXT (one per line) and CSV files.

    Args:
        input_path: Path to the input file.

    Returns:
        List of non-empty query strings.
    """
    input_path = Path(input_path)
    suffix = input_path.suffix.lower()

    if suffix == ".csv":
        queries: list[str] = []
        with open(input_path, encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames:
                # Look for a 'query', 'smiles', 'name', 'cid', or 'compound' column
                target_col: str | None = None
                for col in ["query", "smiles", "name", "cid", "compound", "identifier"]:
                    if col in [c.lower() for c in reader.fieldnames]:
                        for fname in reader.fieldnames:
                            if fname.lower() == col:
                                target_col = fname
                                break
                        break
                if target_col is None:
                    target_col = reader.fieldnames[0]  # Use first column

                for row in reader:
                    val = row.get(target_col, "").strip()
                    if val and not val.startswith("#"):
                        queries.append(val)
            else:
                # No headers, treat as plain text
                f.seek(0)
                for line in f:
                    stripped = line.strip()
                    if stripped and not stripped.startswith("#"):
                        queries.append(stripped)
        return queries
    else:
        # TXT or any other format: one query per line
        queries = []
        with open(input_path, encoding="utf-8") as f:
            for line in f:
                stripped = line.strip()
                if stripped and not stripped.startswith("#"):
                    queries.append(stripped)
        return queries


async def batch_search_async(
    input_path: str | Path,
    sources: list[str] | None = None,
    concurrency: int = 5,
    progress_callback: Callable[[int, int, str], None] | None = None,
) -> tuple[CompoundCollection, list[dict[str, Any]]]:
    """Search for multiple compounds from a file asynchronously.

    Args:
        input_path: Path to the input file (CSV or TXT).
        sources: List of source names to search. Defaults to ["pubchem"].
        concurrency: Maximum concurrent searches.
        progress_callback: Optional callback(completed, total, current_item).

    Returns:
        Tuple of (CompoundCollection with results, list of error dicts).
    """
    from chemfuse.models.collection import CompoundCollection
    from chemfuse.sources import registry

    input_path = Path(input_path)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    queries = _read_queries_from_file(input_path)
    if not queries:
        return CompoundCollection(compounds=[], query="batch", sources=sources or ["pubchem"]), []

    if len(queries) > MAX_BATCH_SIZE:
        raise ValueError(f"Batch size {len(queries)} exceeds maximum {MAX_BATCH_SIZE}")

    active_sources = sources or ["pubchem"]
    total = len(queries)
    results: list[Compound] = []
    errors: list[dict[str, Any]] = []
    semaphore = asyncio.Semaphore(concurrency)
    completed = 0

    async def search_single(query: str) -> None:
        nonlocal completed
        async with semaphore:
            try:
                query_type = _detect_query_type(query)
                for source_name in active_sources:
                    try:
                        adapter = registry.get(source_name)
                        source_results = await adapter.search(query, query_type)
                        results.extend(source_results)
                    except Exception as exc:
                        errors.append({
                            "query": query,
                            "source": source_name,
                            "error": str(exc),
                        })
            except Exception as exc:
                errors.append({"query": query, "source": "unknown", "error": str(exc)})
            finally:
                completed += 1
                if progress_callback:
                    progress_callback(completed, total, query)

    tasks = [search_single(q) for q in queries]
    await asyncio.gather(*tasks)

    # Deduplicate by InChIKey or SMILES to avoid redundant results
    seen: set[str] = set()
    unique: list[Compound] = []
    for c in results:
        key = c.inchikey or c.smiles or str(id(c))
        if key not in seen:
            seen.add(key)
            unique.append(c)

    collection = CompoundCollection(
        compounds=unique,
        query=f"batch:{input_path.name}",
        sources=active_sources,
    )
    return collection, errors


def batch_search(
    input_path: str | Path,
    sources: list[str] | None = None,
    concurrency: int = 5,
    progress_callback: Callable[[int, int, str], None] | None = None,
) -> tuple[CompoundCollection, list[dict[str, Any]]]:
    """Synchronous wrapper for batch_search_async.

    Args:
        input_path: Path to the input file.
        sources: List of source names.
        concurrency: Maximum concurrent searches.
        progress_callback: Optional progress callback.

    Returns:
        Tuple of (CompoundCollection, list of error dicts).
    """
    return asyncio.run(
        batch_search_async(input_path, sources, concurrency, progress_callback)
    )
