"""Tests for shared web utilities (web/_utils.py)."""

from __future__ import annotations

import asyncio

import pandas as pd
import pytest


class TestRunAsync:
    """Tests for _run_async()."""

    def test_runs_simple_coroutine(self) -> None:
        """Executes a simple coroutine and returns its result."""
        from chemfuse.web._utils import _run_async

        async def _coro() -> int:
            return 42

        result = _run_async(_coro())
        assert result == 42

    def test_returns_string_result(self) -> None:
        """Returns string result from a coroutine."""
        from chemfuse.web._utils import _run_async

        async def _coro() -> str:
            return "hello"

        result = _run_async(_coro())
        assert result == "hello"

    def test_returns_none_result(self) -> None:
        """Returns None when coroutine returns None."""
        from chemfuse.web._utils import _run_async

        async def _coro() -> None:
            return None

        result = _run_async(_coro())
        assert result is None

    def test_returns_list_result(self) -> None:
        """Returns list result from a coroutine."""
        from chemfuse.web._utils import _run_async

        async def _coro() -> list:
            return [1, 2, 3]

        result = _run_async(_coro())
        assert result == [1, 2, 3]

    def test_propagates_exception(self) -> None:
        """Propagates exceptions raised by the coroutine."""
        from chemfuse.web._utils import _run_async

        async def _coro() -> None:
            raise ValueError("test error")

        with pytest.raises(ValueError, match="test error"):
            _run_async(_coro())

    def test_runs_coroutine_with_await(self) -> None:
        """Correctly runs a coroutine that itself awaits."""
        from chemfuse.web._utils import _run_async

        async def _inner() -> int:
            await asyncio.sleep(0)
            return 99

        async def _outer() -> int:
            return await _inner()

        result = _run_async(_outer())
        assert result == 99

    def test_runs_when_loop_is_running(self) -> None:
        """Runs coroutine in thread pool when an event loop is already running."""
        from chemfuse.web._utils import _run_async

        async def _test() -> None:
            async def _coro() -> int:
                return 123

            # Call _run_async from within a running event loop
            result = _run_async(_coro())
            assert result == 123

        asyncio.run(_test())


class TestFindSmilesColumn:
    """Tests for _find_smiles_column()."""

    def test_finds_lowercase_smiles(self) -> None:
        """Detects 'smiles' column (lowercase)."""
        from chemfuse.web._utils import _find_smiles_column

        result = _find_smiles_column(["smiles", "name"])
        assert result == "smiles"

    def test_finds_uppercase_smiles(self) -> None:
        """Detects 'SMILES' column (uppercase, case-insensitive)."""
        from chemfuse.web._utils import _find_smiles_column

        result = _find_smiles_column(["SMILES", "name"])
        assert result == "SMILES"

    def test_finds_canonical_smiles(self) -> None:
        """Detects 'canonical_smiles' column."""
        from chemfuse.web._utils import _find_smiles_column

        result = _find_smiles_column(["compound_id", "canonical_smiles"])
        assert result == "canonical_smiles"

    def test_finds_canonical_smiles_uppercase(self) -> None:
        """Detects 'Canonical_SMILES' case-insensitively."""
        from chemfuse.web._utils import _find_smiles_column

        result = _find_smiles_column(["Canonical_SMILES"])
        assert result == "Canonical_SMILES"

    def test_finds_canonicalsmiles_no_underscore(self) -> None:
        """Detects 'canonicalsmiles' (no underscore) column."""
        from chemfuse.web._utils import _find_smiles_column

        result = _find_smiles_column(["canonicalsmiles", "name"])
        assert result == "canonicalsmiles"

    def test_finds_canonicalsmiles_uppercase(self) -> None:
        """Detects 'CanonicalSMILES' (mixed case, no underscore) column."""
        from chemfuse.web._utils import _find_smiles_column

        result = _find_smiles_column(["CanonicalSMILES"])
        assert result == "CanonicalSMILES"

    def test_returns_none_when_no_smiles_column(self) -> None:
        """Returns None when no SMILES-like column is present."""
        from chemfuse.web._utils import _find_smiles_column

        result = _find_smiles_column(["compound_id", "name", "formula"])
        assert result is None

    def test_returns_none_for_empty_list(self) -> None:
        """Returns None for an empty column list."""
        from chemfuse.web._utils import _find_smiles_column

        result = _find_smiles_column([])
        assert result is None

    def test_works_with_pandas_index(self) -> None:
        """Works with a pandas Index object as input."""
        from chemfuse.web._utils import _find_smiles_column

        df = pd.DataFrame({"smiles": ["CCO"], "name": ["ethanol"]})
        result = _find_smiles_column(df.columns)
        assert result == "smiles"

    def test_returns_first_match(self) -> None:
        """Returns the first matching column when multiple matches exist."""
        from chemfuse.web._utils import _find_smiles_column

        result = _find_smiles_column(["smiles", "canonical_smiles"])
        assert result == "smiles"

    def test_does_not_match_partial_names(self) -> None:
        """Does not match column names that only partially contain 'smiles'."""
        from chemfuse.web._utils import _find_smiles_column

        result = _find_smiles_column(["isomeric_smiles", "name"])
        assert result is None
