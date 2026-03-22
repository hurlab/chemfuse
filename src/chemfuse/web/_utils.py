"""Shared utilities for the ChemFuse web UI."""

from __future__ import annotations

import asyncio
import concurrent.futures

import pandas as pd


def _run_async(coro: object) -> object:
    """Run an async coroutine safely from a synchronous context.

    Handles both cases: when no event loop is running (standard case) and
    when an event loop is already running (e.g. inside Jupyter or certain
    async frameworks) by delegating to a thread pool executor.

    Args:
        coro: An awaitable coroutine object to execute.

    Returns:
        The result returned by the coroutine.
    """
    try:
        loop = asyncio.get_running_loop()
    except RuntimeError:
        loop = None

    if loop and loop.is_running():
        with concurrent.futures.ThreadPoolExecutor() as pool:
            return pool.submit(asyncio.run, coro).result()
    return asyncio.run(coro)


def _find_smiles_column(columns: list[str] | pd.Index) -> str | None:
    """Find the SMILES column in a DataFrame by checking common column names.

    Args:
        columns: Column names to search through.

    Returns:
        The matching column name, or None if not found.
    """
    for col in columns:
        if col.lower() in ("smiles", "canonical_smiles", "canonicalsmiles"):
            return col
    return None
