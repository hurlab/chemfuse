"""Tests for the ChemFuse pandas DataFrame accessor (.cf).

Uses real SMILES for aspirin, ibuprofen, and caffeine. No mocks needed for RDKit.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

import chemfuse  # noqa: F401  -- triggers accessor registration
from chemfuse.models.collection import CompoundCollection

# ---------------------------------------------------------------------------
# Test molecules
# ---------------------------------------------------------------------------

ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
IBUPROFEN_SMILES = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
CAFFEINE_SMILES = "Cn1cnc2c1c(=O)n(c(=O)n2C)C"
INVALID_SMILES = "NOT_A_SMILES_STRING"


@pytest.fixture()
def simple_df() -> pd.DataFrame:
    """Three-row DataFrame with a 'smiles' column."""
    return pd.DataFrame(
        {
            "smiles": [ASPIRIN_SMILES, IBUPROFEN_SMILES, CAFFEINE_SMILES],
            "name": ["aspirin", "ibuprofen", "caffeine"],
        }
    )


@pytest.fixture()
def df_with_invalid() -> pd.DataFrame:
    """DataFrame that includes one row with an invalid SMILES string."""
    return pd.DataFrame(
        {
            "smiles": [ASPIRIN_SMILES, INVALID_SMILES, CAFFEINE_SMILES],
        }
    )


@pytest.fixture()
def df_empty() -> pd.DataFrame:
    """Empty DataFrame that still has a 'smiles' column."""
    return pd.DataFrame({"smiles": pd.Series([], dtype=str)})


# ---------------------------------------------------------------------------
# 1. Accessor registration
# ---------------------------------------------------------------------------


class TestAccessorRegistration:
    def test_accessor_registered(self, simple_df: pd.DataFrame) -> None:
        """The .cf accessor is available on DataFrames after chemfuse import."""
        assert hasattr(simple_df, "cf")

    def test_accessor_type(self, simple_df: pd.DataFrame) -> None:
        """The .cf attribute is a ChemFuseAccessor instance."""
        from chemfuse.pandas_ext import ChemFuseAccessor

        assert isinstance(simple_df.cf, ChemFuseAccessor)


# ---------------------------------------------------------------------------
# 2. SMILES column auto-detection
# ---------------------------------------------------------------------------


class TestSmileColDetection:
    @pytest.mark.parametrize(
        "col_name",
        ["smiles", "SMILES", "Smiles", "canonical_smiles", "canonicalsmiles", "smi"],
    )
    def test_detects_various_column_names(self, col_name: str) -> None:
        df = pd.DataFrame({col_name: [ASPIRIN_SMILES]})
        # If detection fails compute_descriptors raises ValueError; success means pass.
        result = df.cf.compute_descriptors()
        assert len(result) == 1

    def test_raises_when_no_smiles_column(self) -> None:
        df = pd.DataFrame({"compound": [ASPIRIN_SMILES], "id": [1]})
        with pytest.raises(ValueError, match="No SMILES column found"):
            df.cf.compute_descriptors()

    def test_raises_on_empty_column_list(self) -> None:
        df = pd.DataFrame()
        with pytest.raises(ValueError, match="No SMILES column found"):
            df.cf.compute_descriptors()


# ---------------------------------------------------------------------------
# 3. compute_descriptors
# ---------------------------------------------------------------------------


class TestComputeDescriptors:
    def test_adds_descriptor_columns(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.compute_descriptors()
        # Should have at least the original columns plus many descriptor columns
        assert len(result.columns) > len(simple_df.columns)
        assert "smiles" in result.columns
        assert "name" in result.columns

    def test_returns_new_dataframe(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.compute_descriptors()
        assert result is not simple_df

    def test_original_not_modified(self, simple_df: pd.DataFrame) -> None:
        original_cols = list(simple_df.columns)
        simple_df.cf.compute_descriptors()
        assert list(simple_df.columns) == original_cols

    def test_row_count_preserved(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.compute_descriptors()
        assert len(result) == len(simple_df)

    def test_known_descriptor_present(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.compute_descriptors()
        # MolWt is always computed by RDKit
        assert "MolWt" in result.columns

    def test_invalid_smiles_produces_nan_row(self, df_with_invalid: pd.DataFrame) -> None:
        result = df_with_invalid.cf.compute_descriptors()
        # Row at index 1 (invalid SMILES) should have NaN for descriptor columns
        descriptor_cols = [c for c in result.columns if c != "smiles"]
        assert descriptor_cols, "Expected descriptor columns to exist"
        invalid_row = result.iloc[1]
        assert all(math.isnan(v) for v in invalid_row[descriptor_cols] if isinstance(v, float))

    def test_empty_dataframe_returns_empty(self, df_empty: pd.DataFrame) -> None:
        result = df_empty.cf.compute_descriptors()
        assert len(result) == 0


# ---------------------------------------------------------------------------
# 4. compute_fingerprints
# ---------------------------------------------------------------------------


class TestComputeFingerprints:
    def test_adds_fingerprint_column(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.compute_fingerprints()
        assert "fingerprint" in result.columns

    def test_fingerprint_is_dict(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.compute_fingerprints()
        fp = result["fingerprint"].iloc[0]
        assert isinstance(fp, dict)
        assert "bits" in fp
        assert "num_bits" in fp

    def test_default_fp_type_is_morgan(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.compute_fingerprints()
        fp = result["fingerprint"].iloc[0]
        assert fp["type"] == "morgan"

    def test_maccs_fp_type(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.compute_fingerprints(fp_type="maccs")
        fp = result["fingerprint"].iloc[0]
        assert fp["num_bits"] == 167

    def test_invalid_smiles_produces_none(self, df_with_invalid: pd.DataFrame) -> None:
        result = df_with_invalid.cf.compute_fingerprints()
        assert result["fingerprint"].iloc[1] is None

    def test_returns_new_dataframe(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.compute_fingerprints()
        assert result is not simple_df


# ---------------------------------------------------------------------------
# 5. filter_druglike
# ---------------------------------------------------------------------------


class TestFilterDruglike:
    def test_lipinski_default(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.filter_druglike()
        assert isinstance(result, pd.DataFrame)
        # Aspirin, ibuprofen, and caffeine all pass Lipinski
        assert len(result) == 3

    @pytest.mark.parametrize("rule", ["lipinski", "veber", "ghose", "egan", "muegge"])
    def test_all_rules_accepted(self, simple_df: pd.DataFrame, rule: str) -> None:
        result = simple_df.cf.filter_druglike(rule=rule)
        assert isinstance(result, pd.DataFrame)

    def test_unknown_rule_raises(self, simple_df: pd.DataFrame) -> None:
        with pytest.raises(ValueError, match="Unknown drug-likeness rule"):
            simple_df.cf.filter_druglike(rule="unknown_rule")

    def test_returns_subset(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.filter_druglike("lipinski")
        assert set(result.index).issubset(set(simple_df.index))

    def test_invalid_smiles_excluded(self, df_with_invalid: pd.DataFrame) -> None:
        result = df_with_invalid.cf.filter_druglike("lipinski")
        # Invalid SMILES row should be excluded
        assert len(result) < len(df_with_invalid)

    def test_original_not_modified(self, simple_df: pd.DataFrame) -> None:
        original_len = len(simple_df)
        simple_df.cf.filter_druglike("lipinski")
        assert len(simple_df) == original_len

    def test_filter_returns_new_dataframe(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.filter_druglike("lipinski")
        assert result is not simple_df


# ---------------------------------------------------------------------------
# 6. predict_admet
# ---------------------------------------------------------------------------


class TestPredictAdmet:
    def test_adds_admet_columns(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.predict_admet()
        admet_cols = [c for c in result.columns if c.startswith("admet_")]
        assert len(admet_cols) > 0

    def test_overall_score_column_present(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.predict_admet()
        assert "admet_overall_score" in result.columns

    def test_overall_score_in_range(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.predict_admet()
        for score in result["admet_overall_score"].dropna():
            assert 0.0 <= score <= 1.0

    def test_row_count_preserved(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.predict_admet()
        assert len(result) == len(simple_df)

    def test_original_columns_preserved(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.predict_admet()
        assert "smiles" in result.columns
        assert "name" in result.columns

    def test_invalid_smiles_produces_nan(self, df_with_invalid: pd.DataFrame) -> None:
        result = df_with_invalid.cf.predict_admet()
        # Invalid row should have NaN for overall score
        score = result["admet_overall_score"].iloc[1]
        assert math.isnan(score)

    def test_returns_new_dataframe(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.predict_admet()
        assert result is not simple_df


# ---------------------------------------------------------------------------
# 7. add_scaffolds
# ---------------------------------------------------------------------------


class TestAddScaffolds:
    def test_adds_scaffold_column(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.add_scaffolds()
        assert "scaffold" in result.columns

    def test_adds_generic_scaffold_column(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.add_scaffolds()
        assert "generic_scaffold" in result.columns

    def test_scaffold_is_string_or_none(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.add_scaffolds()
        for val in result["scaffold"]:
            assert val is None or isinstance(val, str)

    def test_aspirin_scaffold_not_none(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.add_scaffolds()
        # Aspirin has a benzene ring, so it must have a scaffold
        assert result["scaffold"].iloc[0] is not None

    def test_invalid_smiles_scaffold_is_none(self, df_with_invalid: pd.DataFrame) -> None:
        result = df_with_invalid.cf.add_scaffolds()
        assert result["scaffold"].iloc[1] is None

    def test_row_count_preserved(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.add_scaffolds()
        assert len(result) == len(simple_df)

    def test_returns_new_dataframe(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.add_scaffolds()
        assert result is not simple_df


# ---------------------------------------------------------------------------
# 8. standardize
# ---------------------------------------------------------------------------


class TestStandardize:
    def test_adds_standardized_smiles_column(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.standardize()
        assert "standardized_smiles" in result.columns

    def test_standardized_smiles_is_string(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.standardize()
        for val in result["standardized_smiles"].dropna():
            assert isinstance(val, str)
            assert len(val) > 0

    def test_invalid_smiles_produces_none(self, df_with_invalid: pd.DataFrame) -> None:
        result = df_with_invalid.cf.standardize()
        assert result["standardized_smiles"].iloc[1] is None

    def test_returns_new_dataframe(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.standardize()
        assert result is not simple_df

    def test_row_count_preserved(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.standardize()
        assert len(result) == len(simple_df)


# ---------------------------------------------------------------------------
# 9. to_collection
# ---------------------------------------------------------------------------


class TestToCollection:
    def test_returns_compound_collection(self, simple_df: pd.DataFrame) -> None:
        coll = simple_df.cf.to_collection()
        assert isinstance(coll, CompoundCollection)

    def test_collection_length_matches_rows(self, simple_df: pd.DataFrame) -> None:
        coll = simple_df.cf.to_collection()
        assert len(coll) == len(simple_df)

    def test_smiles_mapped_correctly(self, simple_df: pd.DataFrame) -> None:
        coll = simple_df.cf.to_collection()
        assert coll.compounds[0].smiles == ASPIRIN_SMILES
        assert coll.compounds[1].smiles == IBUPROFEN_SMILES

    def test_name_column_mapped(self, simple_df: pd.DataFrame) -> None:
        coll = simple_df.cf.to_collection()
        assert coll.compounds[0].name == "aspirin"

    def test_optional_columns_mapped(self) -> None:
        df = pd.DataFrame(
            {
                "smiles": [ASPIRIN_SMILES],
                "inchikey": ["BSYNRYMUTXBXSQ-UHFFFAOYSA-N"],
                "formula": ["C9H8O4"],
                "cid": [2244],
            }
        )
        coll = df.cf.to_collection()
        compound = coll.compounds[0]
        assert compound.inchikey == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
        assert compound.formula == "C9H8O4"
        assert compound.cid == 2244

    def test_empty_dataframe_returns_empty_collection(self, df_empty: pd.DataFrame) -> None:
        coll = df_empty.cf.to_collection()
        assert len(coll) == 0

    def test_raises_without_smiles_column(self) -> None:
        df = pd.DataFrame({"compound": [ASPIRIN_SMILES]})
        with pytest.raises(ValueError, match="No SMILES column found"):
            df.cf.to_collection()


# ---------------------------------------------------------------------------
# 10. diversity_pick
# ---------------------------------------------------------------------------


class TestDiversityPick:
    def test_returns_n_rows(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.diversity_pick(2)
        assert len(result) == 2

    def test_returns_subset_of_original(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.diversity_pick(2)
        assert set(result.index).issubset(set(simple_df.index))

    def test_pick_all(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.diversity_pick(3)
        assert len(result) == 3

    def test_pick_more_than_available_returns_all_valid(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.diversity_pick(100)
        assert len(result) <= len(simple_df)

    def test_returns_new_dataframe(self, simple_df: pd.DataFrame) -> None:
        result = simple_df.cf.diversity_pick(2)
        assert result is not simple_df


# ---------------------------------------------------------------------------
# 11. tanimoto_matrix
# ---------------------------------------------------------------------------


class TestTanimotoMatrix:
    def test_returns_numpy_array(self, simple_df: pd.DataFrame) -> None:
        mat = simple_df.cf.tanimoto_matrix()
        assert isinstance(mat, np.ndarray)

    def test_matrix_shape(self, simple_df: pd.DataFrame) -> None:
        mat = simple_df.cf.tanimoto_matrix()
        n = len(simple_df)
        assert mat.shape == (n, n)

    def test_diagonal_is_one(self, simple_df: pd.DataFrame) -> None:
        mat = simple_df.cf.tanimoto_matrix()
        for i in range(len(simple_df)):
            assert abs(mat[i, i] - 1.0) < 1e-6

    def test_similarity_in_range(self, simple_df: pd.DataFrame) -> None:
        mat = simple_df.cf.tanimoto_matrix()
        # Ignore NaN values; all others must be in [0, 1]
        valid = mat[~np.isnan(mat)]
        assert np.all(valid >= 0.0)
        assert np.all(valid <= 1.0)

    def test_matrix_is_symmetric(self, simple_df: pd.DataFrame) -> None:
        mat = simple_df.cf.tanimoto_matrix()
        for i in range(len(simple_df)):
            for j in range(len(simple_df)):
                if not (math.isnan(mat[i, j]) or math.isnan(mat[j, i])):
                    assert abs(mat[i, j] - mat[j, i]) < 1e-6

    def test_maccs_fp_type(self, simple_df: pd.DataFrame) -> None:
        mat = simple_df.cf.tanimoto_matrix(fp_type="maccs")
        assert mat.shape == (len(simple_df), len(simple_df))


# ---------------------------------------------------------------------------
# 12. Error handling edge cases
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_nan_smiles_handled_gracefully(self) -> None:
        df = pd.DataFrame({"smiles": [ASPIRIN_SMILES, None, CAFFEINE_SMILES]})
        result = df.cf.compute_descriptors()
        assert len(result) == 3

    def test_empty_string_smiles_handled(self) -> None:
        df = pd.DataFrame({"smiles": [ASPIRIN_SMILES, "", CAFFEINE_SMILES]})
        result = df.cf.compute_descriptors()
        assert len(result) == 3

    def test_invalid_smiles_does_not_raise_in_compute_descriptors(
        self, df_with_invalid: pd.DataFrame
    ) -> None:
        # Should log a warning but not raise
        result = df_with_invalid.cf.compute_descriptors()
        assert len(result) == len(df_with_invalid)

    def test_single_row_dataframe(self) -> None:
        df = pd.DataFrame({"smiles": [ASPIRIN_SMILES]})
        desc = df.cf.compute_descriptors()
        assert len(desc) == 1
        assert "MolWt" in desc.columns

    def test_all_invalid_smiles_filter_returns_empty(self) -> None:
        df = pd.DataFrame({"smiles": [INVALID_SMILES, "ALSO_INVALID"]})
        result = df.cf.filter_druglike("lipinski")
        assert len(result) == 0
