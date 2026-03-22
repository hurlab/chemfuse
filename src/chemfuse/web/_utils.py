"""Shared utilities for the ChemFuse web UI."""

from __future__ import annotations

import asyncio
import concurrent.futures
import hashlib
import json
from typing import Any

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


def _make_params_hash(params: dict[str, Any]) -> str:
    """Return a stable hex digest for a dict of scalar parameters.

    Only JSON-serialisable values are supported; non-serialisable values are
    converted to their string representation before hashing.

    Args:
        params: Mapping of parameter names to their current values.

    Returns:
        A short hexadecimal hash string representing the parameter set.
    """
    serialisable: dict[str, Any] = {}
    for key, value in params.items():
        try:
            json.dumps(value)
            serialisable[key] = value
        except (TypeError, ValueError):
            serialisable[key] = str(value)

    encoded = json.dumps(serialisable, sort_keys=True).encode()
    return hashlib.md5(encoded).hexdigest()  # noqa: S324 – non-security hash


def _params_changed(
    session_state: Any,
    hash_key: str,
    params: dict[str, Any],
) -> bool:
    """Detect whether the current parameter set differs from the stored hash.

    Updates the stored hash in *session_state* whenever a change is detected.

    Args:
        session_state: Streamlit's ``st.session_state`` object.
        hash_key: The key under which the previous hash is stored.
        params: Current parameter values to compare.

    Returns:
        True if the parameters changed since the last call, False otherwise.
    """
    current_hash = _make_params_hash(params)
    previous_hash = session_state.get(hash_key)
    if previous_hash != current_hash:
        session_state[hash_key] = current_hash
        return True
    return False
