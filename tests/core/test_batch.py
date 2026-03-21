"""Tests for chemfuse.core.batch."""

from __future__ import annotations

import csv
from pathlib import Path
from unittest.mock import AsyncMock, patch

import pytest

from chemfuse.core.batch import _detect_query_type, _read_queries_from_file, batch_search
from chemfuse.models.collection import CompoundCollection
from chemfuse.models.compound import Compound, CompoundProperties


@pytest.fixture
def aspirin_compound() -> Compound:
    return Compound(
        cid=2244,
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        name="aspirin",
        sources=["pubchem"],
        properties=CompoundProperties(molecular_weight=180.16),
    )


class TestDetectQueryType:
    def test_integer_string_is_cid(self):
        assert _detect_query_type("2244") == "cid"

    def test_smiles_detected(self):
        assert _detect_query_type("CC(=O)Oc1ccccc1C(=O)O") == "smiles"

    def test_formula_detected(self):
        assert _detect_query_type("C9H8O4") == "formula"

    def test_inchi_detected(self):
        result = _detect_query_type("InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)")
        assert result == "inchi"

    def test_plain_name_is_name(self):
        assert _detect_query_type("aspirin") == "name"
        assert _detect_query_type("acetylsalicylic acid") == "name"

    def test_large_integer_is_cid(self):
        assert _detect_query_type("12345678") == "cid"


class TestReadQueriesFromFile:
    def test_reads_txt_file(self, tmp_path: Path):
        txt_file = tmp_path / "queries.txt"
        txt_file.write_text("aspirin\nibuprofen\ncaffeine\n")
        queries = _read_queries_from_file(str(txt_file))
        assert queries == ["aspirin", "ibuprofen", "caffeine"]

    def test_reads_txt_skips_empty_lines(self, tmp_path: Path):
        txt_file = tmp_path / "queries.txt"
        txt_file.write_text("aspirin\n\nibuprofen\n\n")
        queries = _read_queries_from_file(str(txt_file))
        assert queries == ["aspirin", "ibuprofen"]

    def test_reads_csv_query_column(self, tmp_path: Path):
        csv_file = tmp_path / "queries.csv"
        with open(csv_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["query", "other"])
            writer.writeheader()
            writer.writerow({"query": "aspirin", "other": "x"})
            writer.writerow({"query": "ibuprofen", "other": "y"})
        queries = _read_queries_from_file(str(csv_file))
        assert queries == ["aspirin", "ibuprofen"]

    def test_reads_csv_name_column(self, tmp_path: Path):
        csv_file = tmp_path / "queries.csv"
        with open(csv_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["name", "mw"])
            writer.writeheader()
            writer.writerow({"name": "aspirin", "mw": "180.16"})
            writer.writerow({"name": "ibuprofen", "mw": "206.28"})
        queries = _read_queries_from_file(str(csv_file))
        assert queries == ["aspirin", "ibuprofen"]

    def test_reads_csv_smiles_column(self, tmp_path: Path):
        csv_file = tmp_path / "queries.csv"
        with open(csv_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["smiles"])
            writer.writeheader()
            writer.writerow({"smiles": "CC(=O)Oc1ccccc1C(=O)O"})
        queries = _read_queries_from_file(str(csv_file))
        assert queries == ["CC(=O)Oc1ccccc1C(=O)O"]

    def test_reads_csv_cid_column(self, tmp_path: Path):
        csv_file = tmp_path / "queries.csv"
        with open(csv_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["cid"])
            writer.writeheader()
            writer.writerow({"cid": "2244"})
            writer.writerow({"cid": "3672"})
        queries = _read_queries_from_file(str(csv_file))
        assert queries == ["2244", "3672"]

    def test_raises_for_nonexistent_file(self):
        with pytest.raises((FileNotFoundError, ValueError)):
            _read_queries_from_file("/nonexistent/file.txt")


class TestBatchSearch:
    def test_basic_batch_search(self, tmp_path: Path, aspirin_compound: Compound):
        txt_file = tmp_path / "queries.txt"
        txt_file.write_text("aspirin\n")

        with patch("chemfuse.sources.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(return_value=[aspirin_compound])
            mock_registry.get.return_value = mock_adapter

            collection, errors = batch_search(str(txt_file))

        assert isinstance(collection, CompoundCollection)
        assert isinstance(errors, list)

    def test_returns_errors_on_failure(self, tmp_path: Path):
        txt_file = tmp_path / "queries.txt"
        txt_file.write_text("badcompound\n")

        from chemfuse.core.exceptions import NotFoundError

        with patch("chemfuse.sources.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(
                side_effect=NotFoundError("not found", identifier="badcompound")
            )
            mock_registry.get.return_value = mock_adapter

            collection, errors = batch_search(str(txt_file))

        assert isinstance(errors, list)
        assert len(errors) > 0
        assert "badcompound" in str(errors[0])

    def test_progress_callback_called(self, tmp_path: Path, aspirin_compound: Compound):
        txt_file = tmp_path / "queries.txt"
        txt_file.write_text("aspirin\nibuprofen\n")

        progress_calls = []

        def progress_cb(current: int, total: int, query: str) -> None:
            progress_calls.append((current, total, query))

        with patch("chemfuse.sources.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(return_value=[aspirin_compound])
            mock_registry.get.return_value = mock_adapter

            batch_search(str(txt_file), progress_callback=progress_cb)

        assert len(progress_calls) > 0

    def test_batch_with_custom_sources(self, tmp_path: Path, aspirin_compound: Compound):
        txt_file = tmp_path / "queries.txt"
        txt_file.write_text("aspirin\n")

        with patch("chemfuse.sources.registry") as mock_registry:
            mock_adapter = AsyncMock()
            mock_adapter.search = AsyncMock(return_value=[aspirin_compound])
            mock_registry.get.return_value = mock_adapter

            collection, errors = batch_search(str(txt_file), sources=["pubchem"])

        assert isinstance(collection, CompoundCollection)
