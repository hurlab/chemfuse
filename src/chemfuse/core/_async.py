"""Canonical async-to-sync bridge for ChemFuse."""

from __future__ import annotations

import asyncio


def run_async(coro: object) -> object:
    """Run an async coroutine safely from a synchronous context.

    Handles both cases: when no event loop is running (standard script/CLI
    case) and when a loop is already running (e.g., Jupyter notebooks) by
    delegating to a ThreadPoolExecutor so the inner asyncio.run() gets its
    own fresh event loop.

    Args:
        coro: An awaitable coroutine object to execute.

    Returns:
        The result returned by the coroutine.
    """
    try:
        asyncio.get_running_loop()
    except RuntimeError:
        return asyncio.run(coro)

    # There is an active event loop (e.g. Jupyter). Run the coroutine in a
    # dedicated thread so asyncio.run() can create its own nested loop.
    import concurrent.futures

    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
        future = pool.submit(asyncio.run, coro)
        return future.result()
