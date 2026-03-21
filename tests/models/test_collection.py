"""Tests for chemfuse.models.collection."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from chemfuse.models.collection import CompoundCollection
from chemfuse.models.compound import Compound, CompoundProperties


class TestCompoundCollectionCreation:
    def test_empty_collection(self, empty_collection: CompoundCollection):
        assert len(empty_collection) == 0

    def test_collection_with_compounds(self, sample_collection: CompoundCollection):
        assert len(sample_collection) == 3

    def test_repr(self, sample_collection: CompoundCollection):
        r = repr(sample_collection)
        assert "CompoundCollection" in r
        assert "n=3" in r

    def test_iter(self, sample_collection: CompoundCollection):
        names = [c.name for c in sample_collection]
        assert "aspirin" in names
        assert "ibuprofen" in names

    def test_getitem_integer(self, sample_collection: CompoundCollection):
        first = sample_collection[0]
        assert isinstance(first, Compound)

    def test_getitem_slice(self, sample_collection: CompoundCollection):
        subset = sample_collection[0:2]
        assert isinstance(subset, list)
        assert len(subset) == 2

    def test_len_empty(self, empty_collection: CompoundCollection):
        assert len(empty_collection) == 0

    def test_bool_empty(self, empty_collection: CompoundCollection):
        assert not bool(empty_collection)

    def test_bool_nonempty(self, sample_collection: CompoundCollection):
        assert bool(sample_collection)


class TestFilter:
    def test_filter_by_mw_range(self, sample_collection: CompoundCollection):
        result = sample_collection.filter(mw_range=(175.0, 185.0))
        # Only aspirin (MW=180.16) should match
        names = [c.name for c in result]
        assert "aspirin" in names
        assert "ibuprofen" not in names

    def test_filter_by_logp(self, sample_collection: CompoundCollection):
        result = sample_collection.filter(logp_range=(3.0, 5.0))
        names = [c.name for c in result]
        assert "ibuprofen" in names
        assert "caffeine" not in names  # xlogp=-0.07

    def test_filter_by_tpsa(self, sample_collection: CompoundCollection):
        result = sample_collection.filter(tpsa_range=(0.0, 60.0))
        names = [c.name for c in result]
        assert "ibuprofen" in names  # tpsa=37.3
        assert "aspirin" not in names  # tpsa=63.6

    def test_filter_by_hbd_max(self, sample_collection: CompoundCollection):
        result = sample_collection.filter(hbd_max=0)
        names = [c.name for c in result]
        assert "caffeine" in names  # hbd=0
        assert "aspirin" not in names  # hbd=1

    def test_filter_by_hba_max(self, sample_collection: CompoundCollection):
        result = sample_collection.filter(hba_max=2)
        names = [c.name for c in result]
        assert "ibuprofen" in names  # hba=2

    def test_filter_by_sources(self, sample_collection: CompoundCollection):
        result = sample_collection.filter(sources=["pubchem"])
        assert len(result) == 3  # All from pubchem

    def test_filter_returns_new_collection(self, sample_collection: CompoundCollection):
        result = sample_collection.filter(mw_range=(0, 1000))
        assert result is not sample_collection

    def test_filter_preserves_query(self, sample_collection: CompoundCollection):
        result = sample_collection.filter(mw_range=(0, 1000))
        assert result.query == sample_collection.query

    def test_filter_all_none_returns_all(self, sample_collection: CompoundCollection):
        result = sample_collection.filter()
        assert len(result) == len(sample_collection)

    def test_filter_empty_collection(self, empty_collection: CompoundCollection):
        result = empty_collection.filter(mw_range=(0, 500))
        assert len(result) == 0


class TestFilterDruglike:
    def test_lipinski_rule(self, sample_collection: CompoundCollection):
        result = sample_collection.filter_druglike(rule="lipinski")
        # All test compounds should pass Lipinski
        assert len(result) > 0

    def test_veber_rule(self, sample_collection: CompoundCollection):
        result = sample_collection.filter_druglike(rule="veber")
        assert len(result) > 0

    def test_unknown_rule_returns_all(self, sample_collection: CompoundCollection):
        result = sample_collection.filter_druglike(rule="unknown_rule")
        assert len(result) == len(sample_collection)

    def test_lipinski_filters_violators(self):
        """A compound violating Lipinski should be filtered out."""
        violator = Compound(
            smiles="C" * 100,
            name="big_compound",
            sources=["test"],
            properties=CompoundProperties(
                molecular_weight=600.0,  # > 500 violates Lipinski MW
                xlogp=6.0,  # > 5 violates Lipinski XLogP
                hbd_count=6,  # > 5 violates Lipinski HBD
                hba_count=11,  # > 10 violates Lipinski HBA
            ),
        )
        collection = CompoundCollection(compounds=[violator])
        result = collection.filter_druglike(rule="lipinski")
        assert len(result) == 0


class TestSort:
    def test_sort_by_molecular_weight_ascending(self, sample_collection: CompoundCollection):
        result = sample_collection.sort(by="molecular_weight", ascending=True)
        weights = [c.properties.molecular_weight for c in result if c.properties.molecular_weight]
        assert weights == sorted(weights)

    def test_sort_by_molecular_weight_descending(self, sample_collection: CompoundCollection):
        result = sample_collection.sort(by="molecular_weight", ascending=False)
        weights = [c.properties.molecular_weight for c in result if c.properties.molecular_weight]
        assert weights == sorted(weights, reverse=True)

    def test_sort_returns_new_collection(self, sample_collection: CompoundCollection):
        result = sample_collection.sort()
        assert result is not sample_collection

    def test_sort_by_name(self, sample_collection: CompoundCollection):
        result = sample_collection.sort(by="name", ascending=True)
        names = [c.name for c in result if c.name]
        assert names == sorted(names)

    def test_sort_preserves_none_values(self):
        """Compounds with None sort key should end up at the end."""
        c1 = Compound(smiles="C", name="a", properties=CompoundProperties(molecular_weight=None))
        c2 = Compound(smiles="CC", name="b", properties=CompoundProperties(molecular_weight=30.0))
        collection = CompoundCollection(compounds=[c1, c2])
        result = collection.sort(by="molecular_weight", ascending=True)
        # c2 (30.0) should come before c1 (None)
        assert result[0].name == "b"
        assert result[1].name == "a"


class TestToDataFrame:
    def test_returns_dataframe(self, sample_collection: CompoundCollection):
        df = sample_collection.to_dataframe()
        assert isinstance(df, pd.DataFrame)

    def test_correct_number_of_rows(self, sample_collection: CompoundCollection):
        df = sample_collection.to_dataframe()
        assert len(df) == 3

    def test_empty_collection_returns_empty_df(self, empty_collection: CompoundCollection):
        df = empty_collection.to_dataframe()
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0

    def test_columns_include_smiles(self, sample_collection: CompoundCollection):
        df = sample_collection.to_dataframe()
        assert "smiles" in df.columns

    def test_columns_include_name(self, sample_collection: CompoundCollection):
        df = sample_collection.to_dataframe()
        assert "name" in df.columns

    def test_sources_joined_as_string(self, sample_collection: CompoundCollection):
        df = sample_collection.to_dataframe()
        assert "sources" in df.columns
        assert all(isinstance(v, str) for v in df["sources"])

    def test_no_synonyms_column(self, sample_collection: CompoundCollection):
        df = sample_collection.to_dataframe()
        # synonyms list column should not be in DataFrame
        assert "synonyms" not in df.columns

    def test_all_none_columns_dropped(self):
        c = Compound(smiles="C", sources=["test"])
        collection = CompoundCollection(compounds=[c])
        df = collection.to_dataframe()
        # Columns with all-None values should be dropped
        for col in df.columns:
            assert not df[col].isna().all()


class TestExportMethods:
    def test_to_csv(self, sample_collection: CompoundCollection, tmp_path: Path):
        output = str(tmp_path / "test.csv")
        result = sample_collection.to_csv(output)
        assert Path(result).exists()
        df = pd.read_csv(result)
        assert len(df) == 3

    def test_to_json(self, sample_collection: CompoundCollection, tmp_path: Path):
        import json
        output = str(tmp_path / "test.json")
        result = sample_collection.to_json(output)
        assert Path(result).exists()
        with open(result) as f:
            data = json.load(f)
        assert data is not None

    def test_to_excel_with_openpyxl(self, sample_collection: CompoundCollection, tmp_path: Path):
        output = str(tmp_path / "test.xlsx")
        try:
            result = sample_collection.to_excel(output)
            assert Path(result).exists()
        except ImportError:
            pytest.skip("openpyxl not installed")
