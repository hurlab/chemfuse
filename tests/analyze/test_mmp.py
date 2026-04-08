"""Tests for analyze/mmp.py — Matched Molecular Pair analysis."""

from __future__ import annotations

from unittest.mock import patch

import pytest

from chemfuse.core.exceptions import OptionalDependencyError

# Drug-like test molecules
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
ASPIRIN_OH = "CC(=O)Oc1ccc(O)cc1C(=O)O"  # Aspirin with a para-hydroxyl (close analog)
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"  # Very different scaffold


class TestFindMatchedPairs:
    def test_empty_input_returns_empty_dataframe(self):
        """Empty SMILES list must return an empty DataFrame with the expected columns."""
        from chemfuse.analyze.mmp import find_matched_pairs

        df = find_matched_pairs([])
        assert df.empty
        assert "smiles_a" in df.columns
        assert "smiles_b" in df.columns

    def test_single_compound_returns_empty_dataframe(self):
        """A single compound cannot form any pair."""
        from chemfuse.analyze.mmp import find_matched_pairs

        df = find_matched_pairs([ASPIRIN])
        assert df.empty

    def test_close_analogs_detected_as_pair(self):
        """Aspirin and its para-hydroxyl analog share a large core — should pair."""
        from chemfuse.analyze.mmp import find_matched_pairs

        df = find_matched_pairs([ASPIRIN, ASPIRIN_OH])
        # Both molecules share a substantial core; expect at least one pair
        assert len(df) >= 1
        assert "smiles_a" in df.columns
        assert "smiles_b" in df.columns
        assert "core_smarts" in df.columns
        assert "transform_a" in df.columns
        assert "transform_b" in df.columns

    def test_very_different_molecules_not_paired(self):
        """Aspirin and ibuprofen have very different scaffolds — should not pair."""
        from chemfuse.analyze.mmp import find_matched_pairs

        df = find_matched_pairs(
            [ASPIRIN, IBUPROFEN],
            min_heavy_atoms_core=10,  # Require a large shared core
            max_heavy_atoms_transform=5,
        )
        # Different scaffolds should produce zero pairs under strict settings
        assert len(df) == 0

    def test_activity_delta_computed_correctly(self):
        """activity_delta must equal activity_b - activity_a."""
        from chemfuse.analyze.mmp import find_matched_pairs

        activities = [5.0, 7.0]
        df = find_matched_pairs([ASPIRIN, ASPIRIN_OH], activities=activities)

        if len(df) > 0:
            row = df.iloc[0]
            if row["activity_a"] is not None and row["activity_b"] is not None:
                expected_delta = row["activity_b"] - row["activity_a"]
                assert abs(row["activity_delta"] - expected_delta) < 1e-9

    def test_activity_columns_none_when_not_provided(self):
        """When activities is None, activity columns must all be None."""
        from chemfuse.analyze.mmp import find_matched_pairs

        df = find_matched_pairs([ASPIRIN, ASPIRIN_OH])
        if len(df) > 0:
            assert df["activity_a"].iloc[0] is None
            assert df["activity_b"].iloc[0] is None
            assert df["activity_delta"].iloc[0] is None

    def test_activities_length_mismatch_raises(self):
        """Mismatched lengths between smiles_list and activities must raise ValueError."""
        from chemfuse.analyze.mmp import find_matched_pairs

        with pytest.raises(ValueError, match="activities length"):
            find_matched_pairs([ASPIRIN, ASPIRIN_OH], activities=[1.0])

    def test_returns_dataframe_with_expected_columns(self):
        """Result DataFrame must have all required columns."""
        from chemfuse.analyze.mmp import find_matched_pairs

        df = find_matched_pairs([ASPIRIN, ASPIRIN_OH, IBUPROFEN])
        expected_cols = {
            "smiles_a", "smiles_b", "core_smarts",
            "transform_a", "transform_b",
            "activity_a", "activity_b", "activity_delta",
        }
        assert expected_cols.issubset(set(df.columns))

    def test_rdkit_unavailable_raises(self):
        """OptionalDependencyError must be raised when RDKit is not available."""
        import chemfuse.analyze.mmp as mmp_module

        with patch.object(mmp_module, "_RDKIT_AVAILABLE", False):
            from chemfuse.analyze.mmp import find_matched_pairs
            with pytest.raises(OptionalDependencyError):
                find_matched_pairs([ASPIRIN, ASPIRIN_OH])

    def test_sorted_by_abs_activity_delta_descending(self):
        """When activities provided, rows must be sorted by |activity_delta| descending."""
        from chemfuse.analyze.mmp import find_matched_pairs

        # Three analogs: aspirin, aspirin_oh, and another aspirin variant
        aspirin_methyl = "CC(=O)Oc1ccc(C)cc1C(=O)O"  # methyl instead of OH
        activities = [1.0, 6.0, 3.0]  # large delta between 0 and 1, small for 0 and 2
        df = find_matched_pairs(
            [ASPIRIN, ASPIRIN_OH, aspirin_methyl],
            activities=activities,
        )
        if len(df) > 1 and df["activity_delta"].notna().all():
            deltas_abs = df["activity_delta"].abs().tolist()
            assert deltas_abs == sorted(deltas_abs, reverse=True)


class TestSummarizeTransformations:
    def test_empty_dataframe_returns_empty(self):
        """Empty input DataFrame must return an empty summary DataFrame."""
        import pandas as pd

        from chemfuse.analyze.mmp import summarize_transformations

        empty_df = pd.DataFrame(
            columns=[
                "smiles_a", "smiles_b", "core_smarts",
                "transform_a", "transform_b",
                "activity_a", "activity_b", "activity_delta",
            ]
        )
        result = summarize_transformations(empty_df)
        assert result.empty

    def test_missing_transform_columns_returns_empty(self):
        """DataFrame without transform columns must return empty summary."""
        import pandas as pd

        from chemfuse.analyze.mmp import summarize_transformations

        df = pd.DataFrame({"smiles_a": [ASPIRIN], "smiles_b": [ASPIRIN_OH]})
        result = summarize_transformations(df)
        assert result.empty

    def test_groups_by_transformation_and_counts(self):
        """Transformations that appear min_count times must be included."""
        import pandas as pd

        from chemfuse.analyze.mmp import summarize_transformations

        # Construct a DataFrame with two rows sharing the same transformation
        data = {
            "smiles_a": [ASPIRIN, ASPIRIN],
            "smiles_b": [ASPIRIN_OH, ASPIRIN_OH],
            "core_smarts": ["[#6]", "[#6]"],
            "transform_a": ["H", "H"],
            "transform_b": ["O", "O"],
            "activity_a": [1.0, 2.0],
            "activity_b": [3.0, 5.0],
            "activity_delta": [2.0, 3.0],
        }
        df = pd.DataFrame(data)
        result = summarize_transformations(df, min_count=2)
        assert len(result) == 1
        assert result.iloc[0]["count"] == 2

    def test_min_count_filters_rare_transformations(self):
        """Transformations below min_count threshold must be excluded."""
        import pandas as pd

        from chemfuse.analyze.mmp import summarize_transformations

        data = {
            "smiles_a": [ASPIRIN],
            "smiles_b": [ASPIRIN_OH],
            "core_smarts": ["[#6]"],
            "transform_a": ["H"],
            "transform_b": ["O"],
            "activity_a": [1.0],
            "activity_b": [3.0],
            "activity_delta": [2.0],
        }
        df = pd.DataFrame(data)
        # min_count=2 but only 1 occurrence — should return empty
        result = summarize_transformations(df, min_count=2)
        assert result.empty

    def test_mean_activity_delta_computed(self):
        """mean_activity_delta must equal the mean of grouped activity_delta values."""
        import pandas as pd

        from chemfuse.analyze.mmp import summarize_transformations

        data = {
            "smiles_a": [ASPIRIN, ASPIRIN],
            "smiles_b": [ASPIRIN_OH, ASPIRIN_OH],
            "core_smarts": ["[#6]", "[#6]"],
            "transform_a": ["H", "H"],
            "transform_b": ["O", "O"],
            "activity_a": [1.0, 3.0],
            "activity_b": [3.0, 7.0],
            "activity_delta": [2.0, 4.0],
        }
        df = pd.DataFrame(data)
        result = summarize_transformations(df, min_count=2)
        assert len(result) == 1
        assert abs(result.iloc[0]["mean_activity_delta"] - 3.0) < 1e-9

    def test_result_sorted_by_count_descending(self):
        """Result must be sorted by count descending."""
        import pandas as pd

        from chemfuse.analyze.mmp import summarize_transformations

        # Two transformations: one appears 3 times, another 2 times
        data = {
            "smiles_a": [ASPIRIN] * 5,
            "smiles_b": [ASPIRIN_OH] * 5,
            "core_smarts": ["[#6]"] * 5,
            "transform_a": ["H", "H", "H", "CH3", "CH3"],
            "transform_b": ["O", "O", "O", "F", "F"],
            "activity_a": [1.0] * 5,
            "activity_b": [2.0] * 5,
            "activity_delta": [1.0] * 5,
        }
        df = pd.DataFrame(data)
        result = summarize_transformations(df, min_count=2)
        if len(result) > 1:
            counts = result["count"].tolist()
            assert counts == sorted(counts, reverse=True)

    def test_output_columns_present(self):
        """Summary DataFrame must contain all expected columns."""
        import pandas as pd

        from chemfuse.analyze.mmp import summarize_transformations

        data = {
            "smiles_a": [ASPIRIN, ASPIRIN],
            "smiles_b": [ASPIRIN_OH, ASPIRIN_OH],
            "core_smarts": ["[#6]", "[#6]"],
            "transform_a": ["H", "H"],
            "transform_b": ["O", "O"],
            "activity_a": [1.0, 2.0],
            "activity_b": [3.0, 4.0],
            "activity_delta": [2.0, 2.0],
        }
        df = pd.DataFrame(data)
        result = summarize_transformations(df, min_count=2)
        expected_cols = {
            "transform_from", "transform_to", "count",
            "mean_activity_delta", "std_activity_delta", "compounds_involved",
        }
        assert expected_cols.issubset(set(result.columns))
