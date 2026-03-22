"""Tests for chemfuse.cli._async._run_async helper."""

from __future__ import annotations

import asyncio

import pytest


class TestCliRunAsync:
    """Tests for the CLI _run_async() helper."""

    def test_runs_simple_coroutine(self) -> None:
        """Executes a simple coroutine and returns its result."""
        from chemfuse.cli._async import _run_async

        async def _coro() -> int:
            return 42

        assert _run_async(_coro()) == 42

    def test_returns_string_result(self) -> None:
        """Returns string result from a coroutine."""
        from chemfuse.cli._async import _run_async

        async def _coro() -> str:
            return "hello"

        assert _run_async(_coro()) == "hello"

    def test_returns_none_result(self) -> None:
        """Returns None when coroutine returns None."""
        from chemfuse.cli._async import _run_async

        async def _coro() -> None:
            return None

        assert _run_async(_coro()) is None

    def test_returns_list_result(self) -> None:
        """Returns list result from a coroutine."""
        from chemfuse.cli._async import _run_async

        async def _coro() -> list:
            return [1, 2, 3]

        assert _run_async(_coro()) == [1, 2, 3]

    def test_propagates_exception(self) -> None:
        """Propagates exceptions raised by the coroutine."""
        from chemfuse.cli._async import _run_async

        async def _coro() -> None:
            raise ValueError("test error")

        with pytest.raises(ValueError, match="test error"):
            _run_async(_coro())

    def test_runs_coroutine_with_await(self) -> None:
        """Correctly executes a coroutine that itself awaits."""
        from chemfuse.cli._async import _run_async

        async def _inner() -> int:
            await asyncio.sleep(0)
            return 99

        async def _outer() -> int:
            return await _inner()

        assert _run_async(_outer()) == 99

    def test_runs_when_loop_is_already_running(self) -> None:
        """Runs coroutine via thread pool when an event loop is already running.

        This verifies the Jupyter notebook / embedded event-loop scenario:
        calling _run_async() from inside a running loop must not raise
        'This event loop is already running'.
        """
        from chemfuse.cli._async import _run_async

        async def _test() -> None:
            async def _coro() -> int:
                return 123

            result = _run_async(_coro())
            assert result == 123

        asyncio.run(_test())
