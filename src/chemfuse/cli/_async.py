"""Async execution helper for ChemFuse CLI commands.

Provides a safe wrapper for running async coroutines from synchronous
Click command handlers, including environments where an event loop is
already running (e.g., Jupyter notebooks).
"""

from __future__ import annotations

import asyncio
import concurrent.futures


def _run_async(coro: object) -> object:
    """Run an async coroutine safely from a synchronous context.

    Handles both cases: when no event loop is running (standard CLI case)
    and when an event loop is already running (e.g., inside a Jupyter
    notebook) by delegating to a thread pool executor.

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
