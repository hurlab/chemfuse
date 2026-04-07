"""Async execution helper for ChemFuse CLI commands.

Provides a safe wrapper for running async coroutines from synchronous
Click command handlers, including environments where an event loop is
already running (e.g., Jupyter notebooks).
"""

from __future__ import annotations

from chemfuse.core._async import run_async


def _run_async(coro: object) -> object:
    """Run an async coroutine safely from a synchronous context.

    Delegates to the canonical :func:`chemfuse.core._async.run_async`
    implementation.

    Args:
        coro: An awaitable coroutine object to execute.

    Returns:
        The result returned by the coroutine.
    """
    return run_async(coro)
