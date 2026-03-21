"""Tests for chemfuse.core.export."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from chemfuse.core.export import export_csv, export_excel, export_json, export_sdf
from chemfuse.models.collection import CompoundCollection
from chemfuse.models.compound import Compound, CompoundProperties


@pytest.fixture
def sample_collection() -> CompoundCollection:
    """Small CompoundCollection for export tests."""
    compounds = [
        Compound(
            cid=2244,
            smiles="CC(=O)Oc1ccccc1C(=O)O",
            name="aspirin",
            formula="C9H8O4",
            sources=["pubchem"],
            properties=CompoundProperties(molecular_weight=180.16, xlogp=1.2),
        ),
        Compound(
            cid=3672,
            smiles="CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            name="ibuprofen",
            formula="C13H18O2",
            sources=["pubchem"],
            properties=CompoundProperties(molecular_weight=206.28, xlogp=3.5),
        ),
    ]
    return CompoundCollection(compounds=compounds, query="test", sources=["pubchem"])


@pytest.fixture
def empty_collection() -> CompoundCollection:
    """Empty CompoundCollection."""
    return CompoundCollection()


class TestExportCSV:
    def test_creates_file(self, sample_collection: CompoundCollection, tmp_path: Path):
        output = str(tmp_path / "output.csv")
        result = export_csv(sample_collection, output)
        assert Path(result).exists()

    def test_csv_has_correct_rows(self, sample_collection: CompoundCollection, tmp_path: Path):
        output = str(tmp_path / "output.csv")
        export_csv(sample_collection, output)
        df = pd.read_csv(output)
        assert len(df) == 2

    def test_csv_has_name_column(self, sample_collection: CompoundCollection, tmp_path: Path):
        output = str(tmp_path / "output.csv")
        export_csv(sample_collection, output)
        df = pd.read_csv(output)
        assert "name" in df.columns

    def test_csv_has_smiles_column(self, sample_collection: CompoundCollection, tmp_path: Path):
        output = str(tmp_path / "output.csv")
        export_csv(sample_collection, output)
        df = pd.read_csv(output)
        assert "smiles" in df.columns

    def test_empty_collection_creates_file(self, empty_collection: CompoundCollection, tmp_path: Path):
        output = str(tmp_path / "empty.csv")
        result = export_csv(empty_collection, output)
        assert Path(result).exists()

    def test_returns_path_string(self, sample_collection: CompoundCollection, tmp_path: Path):
        output = str(tmp_path / "output.csv")
        result = export_csv(sample_collection, output)
        assert isinstance(result, (str, Path))

    def test_creates_parent_dirs(self, sample_collection: CompoundCollection, tmp_path: Path):
        output = str(tmp_path / "subdir" / "output.csv")
        export_csv(sample_collection, output)
        assert Path(output).exists()


class TestExportJSON:
    def test_creates_file(self, sample_collection: CompoundCollection, tmp_path: Path):
        output = str(tmp_path / "output.json")
        result = export_json(sample_collection, output)
        assert Path(result).exists()

    def test_json_is_valid(self, sample_collection: CompoundCollection, tmp_path: Path):
        output = str(tmp_path / "output.json")
        export_json(sample_collection, output)
        with open(output) as f:
            data = json.load(f)
        assert isinstance(data, (list, dict))

    def test_json_has_compounds(self, sample_collection: CompoundCollection, tmp_path: Path):
        output = str(tmp_path / "output.json")
        export_json(sample_collection, output)
        with open(output) as f:
            data = json.load(f)
        # Should have 2 entries
        if isinstance(data, list):
            assert len(data) == 2
        elif isinstance(data, dict) and "compounds" in data:
            assert len(data["compounds"]) == 2

    def test_empty_collection(self, empty_collection: CompoundCollection, tmp_path: Path):
        output = str(tmp_path / "empty.json")
        result = export_json(empty_collection, output)
        assert Path(result).exists()

    def test_json_strips_none_values(self, sample_collection: CompoundCollection, tmp_path: Path):
        output = str(tmp_path / "output.json")
        export_json(sample_collection, output)
        with open(output) as f:
            content = f.read()
        assert "null" not in content


class TestRecordsToDataFrame:
    def test_empty_returns_empty_df(self):
        from chemfuse.core.export import records_to_dataframe
        df = records_to_dataframe([])
        assert len(df) == 0

    def test_converts_records(self):
        from chemfuse.core.export import records_to_dataframe
        records = [{"name": "aspirin", "mw": 180.16}, {"name": "ibuprofen", "mw": 206.28}]
        df = records_to_dataframe(records)
        assert len(df) == 2
        assert "name" in df.columns

    def test_drops_all_none_columns(self):
        from chemfuse.core.export import records_to_dataframe
        records = [{"name": "aspirin", "empty": None}, {"name": "ibuprofen", "empty": None}]
        df = records_to_dataframe(records)
        assert "empty" not in df.columns


class TestExportExcel:
    def test_requires_openpyxl(self, sample_collection: CompoundCollection, tmp_path: Path):
        """Test that export_excel works when openpyxl is installed, or raises ImportError."""
        output = str(tmp_path / "output.xlsx")
        try:
            result = export_excel(sample_collection, output)
            assert Path(result).exists()
        except ImportError as e:
            assert "openpyxl" in str(e).lower() or "excel" in str(e).lower()

    def test_import_error_message(self, tmp_path: Path):
        """If openpyxl is missing, error should mention install instructions."""

        # Only test if openpyxl is not installed
        try:
            import openpyxl  # noqa: F401
            pytest.skip("openpyxl is installed")
        except ImportError:
            output = str(tmp_path / "output.xlsx")
            collection = CompoundCollection()
            with pytest.raises(ImportError) as exc_info:
                export_excel(collection, output)
            assert "pip install" in str(exc_info.value).lower() or "chemfuse" in str(exc_info.value).lower()


class TestExportExcelInstalled:
    def test_export_excel_with_openpyxl_creates_file(
        self, sample_collection: CompoundCollection, tmp_path: Path
    ):
        """Test excel export when openpyxl is available."""
        pytest.importorskip("openpyxl")
        output = str(tmp_path / "output.xlsx")
        result = export_excel(sample_collection, output)
        assert Path(result).exists()

    def test_export_excel_custom_sheet_name(
        self, sample_collection: CompoundCollection, tmp_path: Path
    ):
        """Test excel export with custom sheet name."""
        pytest.importorskip("openpyxl")
        output = str(tmp_path / "output.xlsx")
        result = export_excel(sample_collection, output, sheet_name="MySheet")
        assert Path(result).exists()

    def test_export_excel_empty_collection(self, empty_collection: CompoundCollection, tmp_path: Path):
        """Test excel export with empty collection."""
        pytest.importorskip("openpyxl")
        output = str(tmp_path / "empty.xlsx")
        result = export_excel(empty_collection, output)
        assert Path(result).exists()


class TestExportSDF:
    def test_requires_rdkit(self, sample_collection: CompoundCollection, tmp_path: Path):
        """Test that export_sdf works when rdkit is installed, or raises ImportError."""
        output = str(tmp_path / "output.sdf")
        try:
            result = export_sdf(sample_collection, output)
            assert Path(result).exists()
        except ImportError as e:
            assert "rdkit" in str(e).lower()

    def test_import_error_message(self, tmp_path: Path):
        """If rdkit is missing, error should mention install instructions."""
        try:
            from rdkit import Chem  # noqa: F401
            pytest.skip("rdkit is installed")
        except ImportError:
            output = str(tmp_path / "output.sdf")
            collection = CompoundCollection()
            with pytest.raises(ImportError) as exc_info:
                export_sdf(collection, output)
            assert "rdkit" in str(exc_info.value).lower()
