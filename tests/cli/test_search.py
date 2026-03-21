"""Tests for chemfuse.cli.commands.search."""

from __future__ import annotations

from unittest.mock import AsyncMock, patch

import pytest
from click.testing import CliRunner

from chemfuse.cli.main import cli
from chemfuse.core.exceptions import NotFoundError
from chemfuse.models.collection import CompoundCollection
from chemfuse.models.compound import Compound, CompoundProperties


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


@pytest.fixture
def aspirin_collection() -> CompoundCollection:
    return CompoundCollection(
        compounds=[
            Compound(
                cid=2244,
                smiles="CC(=O)Oc1ccccc1C(=O)O",
                name="aspirin",
                formula="C9H8O4",
                sources=["pubchem"],
                properties=CompoundProperties(molecular_weight=180.16),
            )
        ],
        query="aspirin",
        sources=["pubchem"],
    )


@pytest.fixture
def empty_collection() -> CompoundCollection:
    return CompoundCollection()


class TestSearchCommand:
    def test_search_basic(self, runner: CliRunner, aspirin_collection: CompoundCollection):
        with patch("chemfuse.search_async", new_callable=AsyncMock) as mock_search:
            mock_search.return_value = aspirin_collection
            result = runner.invoke(cli, ["search", "aspirin"])
        assert result.exit_code == 0

    def test_search_not_found_exits_1(self, runner: CliRunner):
        with patch("chemfuse.search_async", new_callable=AsyncMock) as mock_search:
            mock_search.side_effect = NotFoundError("not found")
            result = runner.invoke(cli, ["search", "nonexistent_compound_xyz"])
        assert result.exit_code != 0

    def test_search_empty_results_exits_0(self, runner: CliRunner, empty_collection: CompoundCollection):
        with patch("chemfuse.search_async", new_callable=AsyncMock) as mock_search:
            mock_search.return_value = empty_collection
            result = runner.invoke(cli, ["search", "nothing"])
        assert result.exit_code == 0

    def test_search_with_type_option(self, runner: CliRunner, aspirin_collection: CompoundCollection):
        with patch("chemfuse.search_async", new_callable=AsyncMock) as mock_search:
            mock_search.return_value = aspirin_collection
            result = runner.invoke(cli, ["search", "2244", "--type", "cid"])
        assert result.exit_code == 0
        call_kwargs = mock_search.call_args[1]
        assert call_kwargs.get("query_type") == "cid"

    def test_search_with_limit(self, runner: CliRunner, aspirin_collection: CompoundCollection):
        with patch("chemfuse.search_async", new_callable=AsyncMock) as mock_search:
            mock_search.return_value = aspirin_collection
            result = runner.invoke(cli, ["search", "aspirin", "--limit", "5"])
        assert result.exit_code == 0
        call_kwargs = mock_search.call_args[1]
        assert call_kwargs.get("limit") == 5

    def test_search_with_source(self, runner: CliRunner, aspirin_collection: CompoundCollection):
        with patch("chemfuse.search_async", new_callable=AsyncMock) as mock_search:
            mock_search.return_value = aspirin_collection
            result = runner.invoke(cli, ["search", "aspirin", "--source", "pubchem"])
        assert result.exit_code == 0

    def test_search_json_format(self, runner: CliRunner, aspirin_collection: CompoundCollection):
        with patch("chemfuse.search_async", new_callable=AsyncMock) as mock_search:
            mock_search.return_value = aspirin_collection
            result = runner.invoke(cli, ["search", "aspirin", "--format", "json"])
        assert result.exit_code == 0
        assert "{" in result.output or "[" in result.output

    def test_search_csv_format(self, runner: CliRunner, aspirin_collection: CompoundCollection):
        with patch("chemfuse.search_async", new_callable=AsyncMock) as mock_search:
            mock_search.return_value = aspirin_collection
            result = runner.invoke(cli, ["search", "aspirin", "--format", "csv"])
        assert result.exit_code == 0

    def test_search_output_to_file(self, runner: CliRunner, aspirin_collection: CompoundCollection, tmp_path):
        output_file = str(tmp_path / "results.csv")
        with patch("chemfuse.search_async", new_callable=AsyncMock) as mock_search:
            mock_search.return_value = aspirin_collection
            result = runner.invoke(cli, ["search", "aspirin", "--output", output_file])
        assert result.exit_code == 0


class TestCLIVersion:
    def test_version_option(self, runner: CliRunner):
        result = runner.invoke(cli, ["--version"])
        assert result.exit_code == 0
        assert "0.1.0" in result.output or "chemfuse" in result.output.lower()

    def test_help_option(self, runner: CliRunner):
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "search" in result.output

    def test_search_help(self, runner: CliRunner):
        result = runner.invoke(cli, ["search", "--help"])
        assert result.exit_code == 0
        assert "QUERY" in result.output or "query" in result.output.lower()
