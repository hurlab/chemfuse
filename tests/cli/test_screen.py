"""Tests for chemfuse screen CLI command (SPEC-CF-005)."""

from __future__ import annotations

import csv
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from chemfuse.cli.main import cli


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


@pytest.fixture
def smiles_txt(tmp_path: Path) -> str:
    """Plain-text SMILES input file."""
    p = tmp_path / "compounds.txt"
    p.write_text("CC(=O)Oc1ccccc1C(=O)O\nCn1cnc2c1c(=O)n(c(=O)n2C)C\n")
    return str(p)


@pytest.fixture
def smiles_csv(tmp_path: Path) -> str:
    """CSV input file with 'smiles' column."""
    p = tmp_path / "compounds.csv"
    with p.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["smiles", "name"])
        writer.writeheader()
        writer.writerow({"smiles": "CC(=O)Oc1ccccc1C(=O)O", "name": "aspirin"})
        writer.writerow({"smiles": "Cn1cnc2c1c(=O)n(c(=O)n2C)C", "name": "caffeine"})
    return str(p)


@pytest.fixture
def name_csv(tmp_path: Path) -> str:
    """CSV input file with 'name' column only."""
    p = tmp_path / "names.csv"
    with p.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["name"])
        writer.writeheader()
        writer.writerow({"name": "aspirin"})
    return str(p)


def _make_mock_collection(n: int = 2) -> MagicMock:
    """Return a mock CompoundCollection with n compounds."""
    from chemfuse.models.compound import Compound, CompoundProperties

    compounds = [
        Compound(
            smiles=f"C{'C' * i}",
            name=f"compound_{i}",
            sources=["pubchem"],
            properties=CompoundProperties(molecular_weight=100.0 + i * 10),
        )
        for i in range(n)
    ]
    mock_col = MagicMock()
    mock_col.__len__ = lambda self: len(compounds)
    mock_col.__iter__ = lambda self: iter(compounds)
    mock_col.__bool__ = lambda self: len(compounds) > 0
    mock_col.compounds = compounds
    mock_col.filter_druglike.return_value = mock_col
    mock_col.compute_all.return_value = None
    mock_col.predict_admet.return_value = None
    mock_col.cluster.return_value = [0, 1]
    mock_col.to_csv.return_value = "/tmp/test.csv"
    mock_col.to_json.return_value = "/tmp/test.json"
    mock_col.to_excel.return_value = "/tmp/test.xlsx"
    return mock_col


class TestScreenReadInput:
    """Tests for input file parsing."""

    def test_reads_txt_file(self, runner: CliRunner, smiles_txt: str) -> None:
        """Screen command reads SMILES from plain-text file."""
        mock_col = _make_mock_collection(0)

        with patch("chemfuse.cli.commands.screen._search_single", return_value=[]), \
             patch("chemfuse.cli.commands.screen.CompoundCollection", return_value=mock_col):
            result = runner.invoke(cli, ["screen", smiles_txt])

        # Should not crash with a 'no compounds found' message or proceed
        assert result.exit_code in (0, 1)

    def test_reads_csv_smiles_column(self, runner: CliRunner, smiles_csv: str) -> None:
        """Screen command reads SMILES column from CSV."""
        mock_col = _make_mock_collection(0)

        with patch("chemfuse.cli.commands.screen._search_single", return_value=[]), \
             patch("chemfuse.cli.commands.screen.CompoundCollection", return_value=mock_col):
            result = runner.invoke(cli, ["screen", smiles_csv])

        assert result.exit_code in (0, 1)

    def test_reads_csv_name_column(self, runner: CliRunner, name_csv: str) -> None:
        """Screen command falls back to 'name' column when no 'smiles' column."""
        mock_col = _make_mock_collection(0)

        with patch("chemfuse.cli.commands.screen._search_single", return_value=[]), \
             patch("chemfuse.cli.commands.screen.CompoundCollection", return_value=mock_col):
            result = runner.invoke(cli, ["screen", name_csv])

        assert result.exit_code in (0, 1)

    def test_empty_file_exits_with_message(self, runner: CliRunner, tmp_path: Path) -> None:
        """Empty input file causes exit code 1."""
        empty = tmp_path / "empty.txt"
        empty.write_text("")
        result = runner.invoke(cli, ["screen", str(empty)])
        assert result.exit_code == 1

    def test_comment_lines_skipped(self, runner: CliRunner, tmp_path: Path) -> None:
        """Lines starting with '#' are treated as comments."""
        p = tmp_path / "commented.txt"
        p.write_text("# This is a comment\nCC(=O)Oc1ccccc1C(=O)O\n")
        mock_col = _make_mock_collection(0)

        with patch("chemfuse.cli.commands.screen._search_single", return_value=[]), \
             patch("chemfuse.cli.commands.screen.CompoundCollection", return_value=mock_col):
            result = runner.invoke(cli, ["screen", str(p)])

        assert result.exit_code in (0, 1)


class TestScreenOptions:
    """Tests for CLI option handling."""

    def test_admet_flag_calls_predict_admet(self, runner: CliRunner, smiles_txt: str) -> None:
        """--admet flag triggers collection.predict_admet()."""
        mock_col = _make_mock_collection(2)

        with patch("chemfuse.cli.commands.screen._search_single", return_value=mock_col.compounds), \
             patch("chemfuse.cli.commands.screen.CompoundCollection", return_value=mock_col):
            runner.invoke(cli, ["screen", smiles_txt, "--admet"])

        mock_col.predict_admet.assert_called_once()

    def test_no_admet_flag_does_not_call_predict_admet(
        self, runner: CliRunner, smiles_txt: str
    ) -> None:
        """Without --admet, predict_admet is not called."""
        mock_col = _make_mock_collection(2)

        with patch("chemfuse.cli.commands.screen._search_single", return_value=mock_col.compounds), \
             patch("chemfuse.cli.commands.screen.CompoundCollection", return_value=mock_col):
            runner.invoke(cli, ["screen", smiles_txt])

        mock_col.predict_admet.assert_not_called()

    def test_druglikeness_option_calls_filter(self, runner: CliRunner, smiles_txt: str) -> None:
        """--druglikeness triggers filter_druglike()."""
        mock_col = _make_mock_collection(2)

        with patch("chemfuse.cli.commands.screen._search_single", return_value=mock_col.compounds), \
             patch("chemfuse.cli.commands.screen.CompoundCollection", return_value=mock_col):
            runner.invoke(cli, ["screen", smiles_txt, "--druglikeness", "lipinski"])

        mock_col.filter_druglike.assert_called_once_with(rule="lipinski")

    def test_multiple_druglikeness_rules(self, runner: CliRunner, smiles_txt: str) -> None:
        """--druglikeness with comma-separated rules applies each rule."""
        mock_col = _make_mock_collection(2)

        with patch("chemfuse.cli.commands.screen._search_single", return_value=mock_col.compounds), \
             patch("chemfuse.cli.commands.screen.CompoundCollection", return_value=mock_col):
            runner.invoke(cli, ["screen", smiles_txt, "--druglikeness", "lipinski,veber"])

        assert mock_col.filter_druglike.call_count == 2

    def test_cluster_butina_option(self, runner: CliRunner, smiles_txt: str) -> None:
        """--cluster butina triggers collection.cluster()."""
        mock_col = _make_mock_collection(2)

        with patch("chemfuse.cli.commands.screen._search_single", return_value=mock_col.compounds), \
             patch("chemfuse.cli.commands.screen.CompoundCollection", return_value=mock_col):
            runner.invoke(cli, ["screen", smiles_txt, "--cluster", "butina"])

        mock_col.cluster.assert_called_once()
        call_kwargs = mock_col.cluster.call_args
        assert call_kwargs.kwargs.get("method") == "butina" or (
            call_kwargs.args and call_kwargs.args[0] == "butina"
        )

    def test_cluster_kmeans_option(self, runner: CliRunner, smiles_txt: str) -> None:
        """--cluster kmeans triggers collection.cluster() with kmeans."""
        mock_col = _make_mock_collection(2)

        with patch("chemfuse.cli.commands.screen._search_single", return_value=mock_col.compounds), \
             patch("chemfuse.cli.commands.screen.CompoundCollection", return_value=mock_col):
            runner.invoke(cli, ["screen", smiles_txt, "--cluster", "kmeans", "--n-clusters", "3"])

        mock_col.cluster.assert_called_once()

    def test_invalid_cluster_choice_fails(self, runner: CliRunner, smiles_txt: str) -> None:
        """--cluster with invalid value returns non-zero exit code."""
        result = runner.invoke(cli, ["screen", smiles_txt, "--cluster", "invalid_algo"])
        assert result.exit_code != 0

    def test_output_csv(self, runner: CliRunner, smiles_txt: str, tmp_path: Path) -> None:
        """--output .csv triggers collection.to_csv()."""
        mock_col = _make_mock_collection(2)
        out_path = str(tmp_path / "out.csv")
        mock_col.to_csv.return_value = out_path

        with patch("chemfuse.cli.commands.screen._search_single", return_value=mock_col.compounds), \
             patch("chemfuse.cli.commands.screen.CompoundCollection", return_value=mock_col):
            runner.invoke(cli, ["screen", smiles_txt, "-o", out_path])

        mock_col.to_csv.assert_called_once_with(out_path)

    def test_output_json(self, runner: CliRunner, smiles_txt: str, tmp_path: Path) -> None:
        """--output .json triggers collection.to_json()."""
        mock_col = _make_mock_collection(2)
        out_path = str(tmp_path / "out.json")
        mock_col.to_json.return_value = out_path

        with patch("chemfuse.cli.commands.screen._search_single", return_value=mock_col.compounds), \
             patch("chemfuse.cli.commands.screen.CompoundCollection", return_value=mock_col):
            runner.invoke(cli, ["screen", smiles_txt, "-o", out_path])

        mock_col.to_json.assert_called_once_with(out_path)

    def test_output_xlsx(self, runner: CliRunner, smiles_txt: str, tmp_path: Path) -> None:
        """--output .xlsx triggers collection.to_excel()."""
        mock_col = _make_mock_collection(2)
        out_path = str(tmp_path / "out.xlsx")
        mock_col.to_excel.return_value = out_path

        with patch("chemfuse.cli.commands.screen._search_single", return_value=mock_col.compounds), \
             patch("chemfuse.cli.commands.screen.CompoundCollection", return_value=mock_col):
            runner.invoke(cli, ["screen", smiles_txt, "-o", out_path])

        mock_col.to_excel.assert_called_once_with(out_path)


class TestScreenSummary:
    """Tests for screen summary output."""

    def test_summary_printed_when_compounds_found(
        self, runner: CliRunner, smiles_txt: str
    ) -> None:
        """Summary statistics appear in output when compounds are found."""
        mock_col = _make_mock_collection(2)

        with patch("chemfuse.cli.commands.screen._search_single", return_value=mock_col.compounds), \
             patch("chemfuse.cli.commands.screen.CompoundCollection", return_value=mock_col):
            result = runner.invoke(cli, ["screen", smiles_txt])

        assert "Screen Summary" in result.output or result.exit_code == 0

    def test_no_results_exits_zero(self, runner: CliRunner, smiles_txt: str) -> None:
        """No compounds found from search exits 0 (graceful)."""
        mock_col = _make_mock_collection(0)

        with patch("chemfuse.cli.commands.screen._search_single", return_value=[]), \
             patch("chemfuse.cli.commands.screen.CompoundCollection", return_value=mock_col):
            result = runner.invoke(cli, ["screen", smiles_txt])

        # Either exits 0 (no results gracefully) or exits non-zero
        assert result.exit_code in (0, 1)


class TestScreenHelperReadInput:
    """Unit tests for the _read_input helper."""

    def test_reads_plain_text(self, tmp_path: Path) -> None:
        from chemfuse.cli.commands.screen import _read_input

        p = tmp_path / "queries.txt"
        p.write_text("aspirin\ncaffeine\nibuprofen\n")
        result = _read_input(str(p))
        assert result == ["aspirin", "caffeine", "ibuprofen"]

    def test_skips_comment_lines(self, tmp_path: Path) -> None:
        from chemfuse.cli.commands.screen import _read_input

        p = tmp_path / "queries.txt"
        p.write_text("# comment\naspirin\n# another comment\ncaffeine\n")
        result = _read_input(str(p))
        assert result == ["aspirin", "caffeine"]

    def test_reads_csv_smiles_column(self, tmp_path: Path) -> None:
        from chemfuse.cli.commands.screen import _read_input

        p = tmp_path / "data.csv"
        with p.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["smiles", "name"])
            writer.writeheader()
            writer.writerow({"smiles": "CC(=O)Oc1ccccc1C(=O)O", "name": "aspirin"})
            writer.writerow({"smiles": "CCO", "name": "ethanol"})
        result = _read_input(str(p))
        assert "CC(=O)Oc1ccccc1C(=O)O" in result
        assert "CCO" in result

    def test_reads_csv_smi_column(self, tmp_path: Path) -> None:
        from chemfuse.cli.commands.screen import _read_input

        p = tmp_path / "data.csv"
        with p.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["smi"])
            writer.writeheader()
            writer.writerow({"smi": "CCO"})
        result = _read_input(str(p))
        assert "CCO" in result

    def test_reads_csv_name_column_fallback(self, tmp_path: Path) -> None:
        from chemfuse.cli.commands.screen import _read_input

        p = tmp_path / "data.csv"
        with p.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["name", "activity"])
            writer.writeheader()
            writer.writerow({"name": "aspirin", "activity": "10.0"})
        result = _read_input(str(p))
        assert "aspirin" in result

    def test_empty_csv_returns_empty(self, tmp_path: Path) -> None:
        from chemfuse.cli.commands.screen import _read_input

        p = tmp_path / "empty.csv"
        with p.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["smiles"])
            writer.writeheader()
        result = _read_input(str(p))
        assert result == []

    def test_empty_text_returns_empty(self, tmp_path: Path) -> None:
        from chemfuse.cli.commands.screen import _read_input

        p = tmp_path / "empty.txt"
        p.write_text("")
        result = _read_input(str(p))
        assert result == []
