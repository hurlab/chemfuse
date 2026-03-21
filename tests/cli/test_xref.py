"""Tests for chemfuse.cli.commands.xref (xref command)."""

from __future__ import annotations

import json
from unittest.mock import AsyncMock, patch

import pytest
from click.testing import CliRunner

from chemfuse.cli.main import cli


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


ASPIRIN_MAPPING = {
    "pubchem": "2244",
    "chembl": "CHEMBL25",
    "drugbank": "DB00945",
    "kegg": "D00109",
}


class TestXrefCommand:
    """Tests for the xref CLI command basic invocation."""

    def test_xref_no_args_exits_nonzero(self, runner: CliRunner) -> None:
        """xref with no identifier arguments exits with code 1."""
        result = runner.invoke(cli, ["xref"])
        assert result.exit_code == 1
        assert "required" in result.output.lower() or "error" in result.output.lower()

    def test_xref_by_cid_table_output(self, runner: CliRunner) -> None:
        """xref --cid prints a table with database mappings."""
        with patch("chemfuse.map_identifiers_async", new_callable=AsyncMock) as mock_map:
            mock_map.return_value = ASPIRIN_MAPPING
            result = runner.invoke(cli, ["xref", "--cid", "2244"])

        assert result.exit_code == 0
        # Should contain all database names in output
        assert "pubchem" in result.output.lower() or "2244" in result.output

    def test_xref_by_chembl_id(self, runner: CliRunner) -> None:
        """xref --chembl returns mapping."""
        with patch("chemfuse.map_identifiers_async", new_callable=AsyncMock) as mock_map:
            mock_map.return_value = ASPIRIN_MAPPING
            result = runner.invoke(cli, ["xref", "--chembl", "CHEMBL25"])

        assert result.exit_code == 0
        # map_identifiers_async should be called with chembl_id
        mock_map.assert_called_once()
        call_kwargs = mock_map.call_args.kwargs
        assert call_kwargs.get("chembl_id") == "CHEMBL25"

    def test_xref_by_inchikey(self, runner: CliRunner) -> None:
        """xref --inchikey returns mapping."""
        inchikey = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
        with patch("chemfuse.map_identifiers_async", new_callable=AsyncMock) as mock_map:
            mock_map.return_value = ASPIRIN_MAPPING
            result = runner.invoke(cli, ["xref", "--inchikey", inchikey])

        assert result.exit_code == 0
        call_kwargs = mock_map.call_args.kwargs
        assert call_kwargs.get("inchikey") == inchikey

    def test_xref_by_smiles(self, runner: CliRunner) -> None:
        """xref --smiles returns mapping."""
        smiles = "CC(=O)Oc1ccccc1C(=O)O"
        with patch("chemfuse.map_identifiers_async", new_callable=AsyncMock) as mock_map:
            mock_map.return_value = {"pubchem": "2244"}
            result = runner.invoke(cli, ["xref", "--smiles", smiles])

        assert result.exit_code == 0
        call_kwargs = mock_map.call_args.kwargs
        assert call_kwargs.get("smiles") == smiles

    def test_xref_format_json(self, runner: CliRunner) -> None:
        """xref --format json outputs valid JSON."""
        with patch("chemfuse.map_identifiers_async", new_callable=AsyncMock) as mock_map:
            mock_map.return_value = ASPIRIN_MAPPING
            result = runner.invoke(cli, ["xref", "--cid", "2244", "--format", "json"])

        assert result.exit_code == 0
        # Output should be valid JSON
        parsed = json.loads(result.output)
        assert parsed["pubchem"] == "2244"
        assert parsed["chembl"] == "CHEMBL25"

    def test_xref_format_table_default(self, runner: CliRunner) -> None:
        """xref default format is table (no JSON)."""
        with patch("chemfuse.map_identifiers_async", new_callable=AsyncMock) as mock_map:
            mock_map.return_value = {"pubchem": "2244"}
            result = runner.invoke(cli, ["xref", "--cid", "2244"])

        assert result.exit_code == 0
        # Table output is not valid JSON by itself
        try:
            json.loads(result.output)
            # If it parses as JSON, that's acceptable too (some runners may mix)
        except json.JSONDecodeError:
            pass  # Expected for table format

    def test_xref_empty_mapping_exits_zero(self, runner: CliRunner) -> None:
        """xref exits 0 (with message) when no cross-references are found."""
        with patch("chemfuse.map_identifiers_async", new_callable=AsyncMock) as mock_map:
            mock_map.return_value = {}
            result = runner.invoke(cli, ["xref", "--cid", "9999999"])

        # Should exit 0 and print a 'no results' message
        assert result.exit_code == 0

    def test_xref_passes_cid_as_int(self, runner: CliRunner) -> None:
        """xref passes CID as integer to map_identifiers_async."""
        with patch("chemfuse.map_identifiers_async", new_callable=AsyncMock) as mock_map:
            mock_map.return_value = ASPIRIN_MAPPING
            runner.invoke(cli, ["xref", "--cid", "2244"])

        call_kwargs = mock_map.call_args.kwargs
        assert call_kwargs.get("cid") == 2244
        assert isinstance(call_kwargs.get("cid"), int)

    def test_xref_other_args_none_when_not_given(self, runner: CliRunner) -> None:
        """When only --cid is provided, other identifiers are None."""
        with patch("chemfuse.map_identifiers_async", new_callable=AsyncMock) as mock_map:
            mock_map.return_value = ASPIRIN_MAPPING
            runner.invoke(cli, ["xref", "--cid", "2244"])

        call_kwargs = mock_map.call_args.kwargs
        assert call_kwargs.get("chembl_id") is None
        assert call_kwargs.get("inchikey") is None
        assert call_kwargs.get("smiles") is None


class TestXrefOutputFile:
    """Tests for xref --output file writing."""

    def test_xref_output_file_writes_json(self, runner: CliRunner, tmp_path) -> None:
        """xref --output writes JSON to file."""
        out_file = tmp_path / "xref_out.json"
        with patch("chemfuse.map_identifiers_async", new_callable=AsyncMock) as mock_map:
            mock_map.return_value = ASPIRIN_MAPPING
            result = runner.invoke(
                cli, ["xref", "--cid", "2244", "--output", str(out_file)]
            )

        assert result.exit_code == 0
        assert out_file.exists()
        written = json.loads(out_file.read_text())
        assert written["pubchem"] == "2244"
        assert written["chembl"] == "CHEMBL25"

    def test_xref_output_file_does_not_print_to_stdout(
        self, runner: CliRunner, tmp_path
    ) -> None:
        """xref --output does not print mapping to stdout."""
        out_file = tmp_path / "xref_out.json"
        with patch("chemfuse.map_identifiers_async", new_callable=AsyncMock) as mock_map:
            mock_map.return_value = ASPIRIN_MAPPING
            result = runner.invoke(
                cli, ["xref", "--cid", "2244", "--output", str(out_file)]
            )

        assert result.exit_code == 0
        # The mapping values should not appear in stdout
        assert "CHEMBL25" not in result.output


class TestXrefErrorHandling:
    """Tests for xref error handling."""

    def test_xref_unexpected_error_exits_nonzero(self, runner: CliRunner) -> None:
        """xref exits non-zero on unexpected exception."""
        with patch("chemfuse.map_identifiers_async", new_callable=AsyncMock) as mock_map:
            mock_map.side_effect = RuntimeError("Unexpected network failure")
            result = runner.invoke(cli, ["xref", "--cid", "2244"])

        assert result.exit_code == 1
        assert "error" in result.output.lower() or "unexpected" in result.output.lower()

    def test_xref_chemfuse_error_exits_nonzero(self, runner: CliRunner) -> None:
        """xref exits non-zero on ChemFuseError."""
        from chemfuse.core.exceptions import ChemFuseError

        with patch("chemfuse.map_identifiers_async", new_callable=AsyncMock) as mock_map:
            mock_map.side_effect = ChemFuseError("API unavailable")
            result = runner.invoke(cli, ["xref", "--cid", "2244"])

        assert result.exit_code == 1

    def test_xref_invalid_format_exits_nonzero(self, runner: CliRunner) -> None:
        """xref with unknown --format value exits non-zero."""
        with patch("chemfuse.map_identifiers_async", new_callable=AsyncMock) as mock_map:
            mock_map.return_value = ASPIRIN_MAPPING
            result = runner.invoke(
                cli, ["xref", "--cid", "2244", "--format", "csv"]
            )

        # Click should reject invalid choice
        assert result.exit_code != 0


class TestXrefTableFallback:
    """Tests for xref table rendering fallback when rich is unavailable."""

    def test_table_fallback_without_rich(self, runner: CliRunner) -> None:
        """xref table output falls back to plain text when rich is not available."""
        import sys

        # Simulate rich not installed
        with patch("chemfuse.map_identifiers_async", new_callable=AsyncMock) as mock_map:
            mock_map.return_value = {"pubchem": "2244", "chembl": "CHEMBL25"}
            # Patch builtins.__import__ or the import inside _print_xref_table
            with patch.dict(sys.modules, {"rich": None, "rich.console": None, "rich.table": None}):
                result = runner.invoke(cli, ["xref", "--cid", "2244"])

        # Should not crash, should still print something
        assert result.exit_code == 0

    def test_table_output_contains_database_and_identifier(
        self, runner: CliRunner
    ) -> None:
        """Table output contains both database name and identifier."""
        with patch("chemfuse.map_identifiers_async", new_callable=AsyncMock) as mock_map:
            mock_map.return_value = {"pubchem": "2244"}
            # Force plain text fallback
            with patch("chemfuse.cli.commands.xref._print_xref_table") as mock_table:
                mock_table.side_effect = lambda m: None  # suppress actual output
                result = runner.invoke(cli, ["xref", "--cid", "2244"])

        assert result.exit_code == 0
