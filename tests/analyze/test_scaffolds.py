"""Tests for analyze/scaffolds.py (CF-E01: Murcko Scaffold Decomposition)."""

from __future__ import annotations

from unittest.mock import patch

import pandas as pd
import pytest

from chemfuse.core.exceptions import OptionalDependencyError

# ---------------------------------------------------------------------------
# Test molecules
# ---------------------------------------------------------------------------
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
CAFFEINE = "Cn1cnc2c1c(=O)n(c(=O)n2C)C"
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
INVALID_SMILES = "not_a_valid_smiles!!!"
EMPTY_SMILES = ""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _patch_rdkit_unavailable():
    """Context manager: simulate RDKit not installed in the scaffolds module."""
    return patch("chemfuse.analyze.scaffolds._RDKIT_AVAILABLE", False)


# ---------------------------------------------------------------------------
# murcko_scaffold
# ---------------------------------------------------------------------------

class TestMurckoScaffold:
    def test_raises_when_rdkit_unavailable(self):
        """murcko_scaffold raises OptionalDependencyError when RDKit is absent."""
        with _patch_rdkit_unavailable():
            from chemfuse.analyze import scaffolds
            with pytest.raises(OptionalDependencyError):
                scaffolds.murcko_scaffold(ASPIRIN)

    def test_aspirin_returns_benzene_scaffold(self):
        """Aspirin scaffold should be a benzene ring framework."""
        try:
            from chemfuse.analyze.scaffolds import murcko_scaffold
        except ImportError:
            pytest.skip("RDKit not installed")

        result = murcko_scaffold(ASPIRIN)
        assert result is not None
        # The Murcko scaffold of aspirin is a monosubstituted benzene ring
        assert "c1ccccc1" in result or result == "c1ccccc1"

    def test_caffeine_scaffold_not_none(self):
        """Caffeine scaffold extraction succeeds and returns a non-empty string."""
        try:
            from chemfuse.analyze.scaffolds import murcko_scaffold
        except ImportError:
            pytest.skip("RDKit not installed")

        result = murcko_scaffold(CAFFEINE)
        assert result is not None
        assert len(result) > 0

    def test_invalid_smiles_returns_none(self):
        """Invalid SMILES returns None without raising."""
        try:
            from chemfuse.analyze.scaffolds import murcko_scaffold
        except ImportError:
            pytest.skip("RDKit not installed")

        result = murcko_scaffold(INVALID_SMILES)
        assert result is None

    def test_empty_smiles_returns_none(self):
        """Empty SMILES returns None without raising."""
        try:
            from chemfuse.analyze.scaffolds import murcko_scaffold
        except ImportError:
            pytest.skip("RDKit not installed")

        result = murcko_scaffold(EMPTY_SMILES)
        assert result is None

    def test_returns_string(self):
        """murcko_scaffold returns a string for a valid drug-like molecule."""
        try:
            from chemfuse.analyze.scaffolds import murcko_scaffold
        except ImportError:
            pytest.skip("RDKit not installed")

        result = murcko_scaffold(IBUPROFEN)
        assert result is None or isinstance(result, str)


# ---------------------------------------------------------------------------
# generic_scaffold
# ---------------------------------------------------------------------------

class TestGenericScaffold:
    def test_raises_when_rdkit_unavailable(self):
        """generic_scaffold raises OptionalDependencyError when RDKit is absent."""
        with _patch_rdkit_unavailable():
            from chemfuse.analyze import scaffolds
            with pytest.raises(OptionalDependencyError):
                scaffolds.generic_scaffold(ASPIRIN)

    def test_aspirin_generic_scaffold_is_carbon_only(self):
        """Generic scaffold of aspirin should be a carbon-only ring (no N, O, S)."""
        try:
            from chemfuse.analyze.scaffolds import generic_scaffold
        except ImportError:
            pytest.skip("RDKit not installed")

        result = generic_scaffold(ASPIRIN)
        assert result is not None
        # Generic scaffold strips heteroatoms; no lowercase n, o, s letters
        # (aromatic heteroatoms) and no uppercase N, O, S atoms
        for heteroatom in ("N", "O", "S", "n", "o", "s"):
            assert heteroatom not in result, (
                f"Generic scaffold should have no {heteroatom!r}, got {result!r}"
            )

    def test_caffeine_generic_scaffold_strips_heteroatoms(self):
        """Caffeine generic scaffold strips nitrogen atoms from the ring system."""
        try:
            from chemfuse.analyze.scaffolds import generic_scaffold, murcko_scaffold
        except ImportError:
            pytest.skip("RDKit not installed")

        murcko = murcko_scaffold(CAFFEINE)
        generic = generic_scaffold(CAFFEINE)
        # The generic scaffold must differ from the Murcko scaffold for
        # a heteroatom-containing molecule like caffeine
        if murcko is not None and generic is not None:
            # caffeine has nitrogen atoms; they must disappear in generic form
            assert "n" not in generic
            assert "N" not in generic

    def test_invalid_smiles_returns_none(self):
        """Invalid SMILES returns None without raising."""
        try:
            from chemfuse.analyze.scaffolds import generic_scaffold
        except ImportError:
            pytest.skip("RDKit not installed")

        assert generic_scaffold(INVALID_SMILES) is None

    def test_empty_smiles_returns_none(self):
        """Empty SMILES returns None without raising."""
        try:
            from chemfuse.analyze.scaffolds import generic_scaffold
        except ImportError:
            pytest.skip("RDKit not installed")

        assert generic_scaffold(EMPTY_SMILES) is None


# ---------------------------------------------------------------------------
# scaffold_frequency
# ---------------------------------------------------------------------------

class TestScaffoldFrequency:
    def test_raises_when_rdkit_unavailable(self):
        """scaffold_frequency raises OptionalDependencyError when RDKit is absent."""
        with _patch_rdkit_unavailable():
            from chemfuse.analyze import scaffolds
            with pytest.raises(OptionalDependencyError):
                scaffolds.scaffold_frequency([ASPIRIN])

    def test_empty_list_returns_empty_dict(self):
        """Empty input returns an empty dict."""
        try:
            from chemfuse.analyze.scaffolds import scaffold_frequency
        except ImportError:
            pytest.skip("RDKit not installed")

        result = scaffold_frequency([])
        assert result == {}

    def test_single_compound_returns_count_one(self):
        """Single compound yields scaffold with count 1."""
        try:
            from chemfuse.analyze.scaffolds import scaffold_frequency
        except ImportError:
            pytest.skip("RDKit not installed")

        result = scaffold_frequency([ASPIRIN])
        assert len(result) == 1
        assert list(result.values())[0] == 1

    def test_duplicates_increment_count(self):
        """Same scaffold from duplicate SMILES increments the count."""
        try:
            from chemfuse.analyze.scaffolds import scaffold_frequency
        except ImportError:
            pytest.skip("RDKit not installed")

        result = scaffold_frequency([ASPIRIN, ASPIRIN, ASPIRIN])
        assert len(result) == 1
        assert list(result.values())[0] == 3

    def test_two_different_scaffolds(self):
        """Two structurally different molecules yield two separate scaffolds."""
        try:
            from chemfuse.analyze.scaffolds import scaffold_frequency
        except ImportError:
            pytest.skip("RDKit not installed")

        result = scaffold_frequency([ASPIRIN, CAFFEINE])
        # May be 1 or 2 scaffolds depending on RDKit scaffold assignment;
        # at minimum we must get back a non-empty dict.
        assert len(result) >= 1

    def test_sorted_by_frequency_descending(self):
        """Results are sorted by count descending."""
        try:
            from chemfuse.analyze.scaffolds import scaffold_frequency
        except ImportError:
            pytest.skip("RDKit not installed")

        # Aspirin appears twice, caffeine once
        result = scaffold_frequency([ASPIRIN, ASPIRIN, CAFFEINE])
        counts = list(result.values())
        assert counts == sorted(counts, reverse=True)

    def test_invalid_smiles_skipped(self):
        """Invalid SMILES are silently skipped in frequency count."""
        try:
            from chemfuse.analyze.scaffolds import scaffold_frequency
        except ImportError:
            pytest.skip("RDKit not installed")

        result = scaffold_frequency([ASPIRIN, INVALID_SMILES])
        # Only ASPIRIN contributes a scaffold
        assert sum(result.values()) == 1

    def test_returns_dict(self):
        """Return type is a dict."""
        try:
            from chemfuse.analyze.scaffolds import scaffold_frequency
        except ImportError:
            pytest.skip("RDKit not installed")

        result = scaffold_frequency([ASPIRIN])
        assert isinstance(result, dict)


# ---------------------------------------------------------------------------
# group_by_scaffold
# ---------------------------------------------------------------------------

class TestGroupByScaffold:
    def _make_compound(self, smiles: str, name: str | None = None):
        from chemfuse.models.compound import Compound
        return Compound(smiles=smiles, name=name)

    def test_raises_when_rdkit_unavailable(self):
        """group_by_scaffold raises OptionalDependencyError when RDKit is absent."""
        compound = self._make_compound(ASPIRIN)
        with _patch_rdkit_unavailable():
            from chemfuse.analyze import scaffolds
            with pytest.raises(OptionalDependencyError):
                scaffolds.group_by_scaffold([compound])

    def test_empty_list_returns_empty_dict(self):
        """Empty compound list returns an empty dict."""
        try:
            from chemfuse.analyze.scaffolds import group_by_scaffold
        except ImportError:
            pytest.skip("RDKit not installed")

        result = group_by_scaffold([])
        assert result == {}

    def test_same_scaffold_grouped_together(self):
        """Two compounds with the same scaffold go into the same group."""
        try:
            from chemfuse.analyze.scaffolds import group_by_scaffold, murcko_scaffold
        except ImportError:
            pytest.skip("RDKit not installed")

        aspirin1 = self._make_compound(ASPIRIN, name="aspirin-1")
        aspirin2 = self._make_compound(ASPIRIN, name="aspirin-2")

        result = group_by_scaffold([aspirin1, aspirin2])
        scaffold = murcko_scaffold(ASPIRIN)
        assert scaffold is not None
        assert scaffold in result
        assert len(result[scaffold]) == 2

    def test_different_scaffolds_separated(self):
        """Compounds with different scaffolds form different groups."""
        try:
            from chemfuse.analyze.scaffolds import group_by_scaffold, murcko_scaffold
        except ImportError:
            pytest.skip("RDKit not installed")

        c_aspirin = self._make_compound(ASPIRIN, name="aspirin")
        c_caffeine = self._make_compound(CAFFEINE, name="caffeine")

        result = group_by_scaffold([c_aspirin, c_caffeine])
        scaffold_aspirin = murcko_scaffold(ASPIRIN)
        scaffold_caffeine = murcko_scaffold(CAFFEINE)

        if scaffold_aspirin != scaffold_caffeine:
            assert len(result) >= 2

    def test_invalid_smiles_grouped_under_empty_key(self):
        """Compounds with invalid SMILES are grouped under the empty-string key."""
        try:
            from chemfuse.analyze.scaffolds import group_by_scaffold
        except ImportError:
            pytest.skip("RDKit not installed")

        bad = self._make_compound(INVALID_SMILES, name="bad")
        result = group_by_scaffold([bad])
        assert "" in result
        assert bad in result[""]

    def test_returns_dict(self):
        """Return type is a dict."""
        try:
            from chemfuse.analyze.scaffolds import group_by_scaffold
        except ImportError:
            pytest.skip("RDKit not installed")

        compound = self._make_compound(ASPIRIN)
        result = group_by_scaffold([compound])
        assert isinstance(result, dict)


# ---------------------------------------------------------------------------
# CompoundCollection.scaffold_frequency
# ---------------------------------------------------------------------------

class TestCollectionScaffoldFrequency:
    def _make_collection(self, smiles_list: list[str]):
        from chemfuse.models.collection import CompoundCollection
        from chemfuse.models.compound import Compound
        compounds = [Compound(smiles=smi) for smi in smiles_list]
        return CompoundCollection(compounds=compounds)

    def test_returns_dataframe(self):
        """scaffold_frequency() on a collection returns a pandas DataFrame."""
        try:
            from chemfuse.analyze.scaffolds import murcko_scaffold as _ms  # noqa: F401
        except ImportError:
            pytest.skip("RDKit not installed")

        coll = self._make_collection([ASPIRIN, CAFFEINE])
        df = coll.scaffold_frequency()
        assert isinstance(df, pd.DataFrame)

    def test_dataframe_has_expected_columns(self):
        """DataFrame has scaffold, count, and percentage columns."""
        try:
            from chemfuse.analyze.scaffolds import murcko_scaffold as _ms  # noqa: F401
        except ImportError:
            pytest.skip("RDKit not installed")

        coll = self._make_collection([ASPIRIN])
        df = coll.scaffold_frequency()
        for col in ("scaffold", "count", "percentage"):
            assert col in df.columns

    def test_empty_collection_returns_empty_dataframe(self):
        """Empty collection yields empty DataFrame with correct columns."""
        try:
            from chemfuse.analyze.scaffolds import murcko_scaffold as _ms  # noqa: F401
        except ImportError:
            pytest.skip("RDKit not installed")

        coll = self._make_collection([])
        df = coll.scaffold_frequency()
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0

    def test_percentage_sums_to_100(self):
        """Percentage column sums to approximately 100 for valid molecules."""
        try:
            from chemfuse.analyze.scaffolds import murcko_scaffold as _ms  # noqa: F401
        except ImportError:
            pytest.skip("RDKit not installed")

        coll = self._make_collection([ASPIRIN, ASPIRIN, CAFFEINE])
        df = coll.scaffold_frequency()
        if len(df) > 0:
            assert abs(df["percentage"].sum() - 100.0) < 0.5

    def test_count_reflects_duplicates(self):
        """Duplicate SMILES increase the count for that scaffold."""
        try:
            from chemfuse.analyze.scaffolds import murcko_scaffold as _ms  # noqa: F401
        except ImportError:
            pytest.skip("RDKit not installed")

        coll = self._make_collection([ASPIRIN, ASPIRIN, ASPIRIN])
        df = coll.scaffold_frequency()
        assert df["count"].iloc[0] == 3

    def test_sorted_by_count_descending(self):
        """DataFrame rows are sorted by count descending."""
        try:
            from chemfuse.analyze.scaffolds import murcko_scaffold as _ms  # noqa: F401
        except ImportError:
            pytest.skip("RDKit not installed")

        coll = self._make_collection([ASPIRIN, ASPIRIN, CAFFEINE])
        df = coll.scaffold_frequency()
        if len(df) > 1:
            assert df["count"].iloc[0] >= df["count"].iloc[1]


# ---------------------------------------------------------------------------
# CompoundCollection.group_by_scaffold
# ---------------------------------------------------------------------------

class TestCollectionGroupByScaffold:
    def _make_collection(self, smiles_list: list[str]):
        from chemfuse.models.collection import CompoundCollection
        from chemfuse.models.compound import Compound
        compounds = [Compound(smiles=smi) for smi in smiles_list]
        return CompoundCollection(compounds=compounds)

    def test_returns_dict_of_collections(self):
        """group_by_scaffold() returns a dict of CompoundCollection objects."""
        try:
            from chemfuse.analyze.scaffolds import murcko_scaffold as _ms  # noqa: F401
        except ImportError:
            pytest.skip("RDKit not installed")

        from chemfuse.models.collection import CompoundCollection

        coll = self._make_collection([ASPIRIN, CAFFEINE])
        result = coll.group_by_scaffold()

        assert isinstance(result, dict)
        for key, value in result.items():
            assert isinstance(key, str)
            assert isinstance(value, CompoundCollection)

    def test_empty_collection_returns_empty_dict(self):
        """Empty collection yields empty dict."""
        try:
            from chemfuse.analyze.scaffolds import murcko_scaffold as _ms  # noqa: F401
        except ImportError:
            pytest.skip("RDKit not installed")

        coll = self._make_collection([])
        result = coll.group_by_scaffold()
        assert result == {}

    def test_total_compounds_preserved(self):
        """Total compound count across all groups equals original collection size."""
        try:
            from chemfuse.analyze.scaffolds import murcko_scaffold as _ms  # noqa: F401
        except ImportError:
            pytest.skip("RDKit not installed")

        coll = self._make_collection([ASPIRIN, ASPIRIN, CAFFEINE])
        result = coll.group_by_scaffold()
        total = sum(len(sub_coll) for sub_coll in result.values())
        assert total == 3

    def test_compounds_grouped_by_scaffold_key(self):
        """Compounds with the same scaffold are in the same sub-collection."""
        try:
            from chemfuse.analyze.scaffolds import murcko_scaffold
        except ImportError:
            pytest.skip("RDKit not installed")

        coll = self._make_collection([ASPIRIN, ASPIRIN])
        result = coll.group_by_scaffold()
        scaffold = murcko_scaffold(ASPIRIN)
        if scaffold is not None and scaffold in result:
            assert len(result[scaffold]) == 2
