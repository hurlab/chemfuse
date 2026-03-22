"""Tests for R-group decomposition and SAR table generation (CF-E07)."""

from __future__ import annotations

import pytest

# Skip entire module if RDKit is not available
pytest.importorskip("rdkit", reason="RDKit is required for R-group tests")

import pandas as pd

from chemfuse.analyze.rgroup import decompose_rgroups, rgroup_sar_table
from chemfuse.models.bioactivity import Bioactivity
from chemfuse.models.collection import CompoundCollection
from chemfuse.models.compound import Compound

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

BENZENE_CORE = "c1ccccc1"  # benzene core SMARTS
TOLUENE = "Cc1ccccc1"      # methyl substituent
ANILINE = "Nc1ccccc1"      # amino substituent
PHENOL = "Oc1ccccc1"       # hydroxyl substituent
CYCLOHEXANE = "C1CCCCC1"   # does not match benzene core


# ---------------------------------------------------------------------------
# decompose_rgroups tests
# ---------------------------------------------------------------------------


class TestDecomposeRgroups:
    """Tests for the decompose_rgroups function."""

    def test_three_substituted_benzenes_returns_r1(self) -> None:
        """Three monosubstituted benzenes should produce a DataFrame with R1."""
        smiles = [TOLUENE, ANILINE, PHENOL]
        df = decompose_rgroups(smiles, BENZENE_CORE)

        assert not df.empty
        assert "smiles" in df.columns
        assert "R1" in df.columns
        assert len(df) == 3

    def test_smiles_column_matches_input(self) -> None:
        """The smiles column should contain the original SMILES strings."""
        smiles = [TOLUENE, ANILINE]
        df = decompose_rgroups(smiles, BENZENE_CORE)

        assert set(df["smiles"]) == {TOLUENE, ANILINE}

    def test_unmatched_compound_excluded(self) -> None:
        """Cyclohexane does not match the benzene core and must be excluded."""
        smiles = [TOLUENE, CYCLOHEXANE, PHENOL]
        df = decompose_rgroups(smiles, BENZENE_CORE)

        assert CYCLOHEXANE not in df["smiles"].tolist()
        assert len(df) == 2

    def test_invalid_core_smarts_raises_value_error(self) -> None:
        """An invalid SMARTS string should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid core SMARTS"):
            decompose_rgroups([TOLUENE], "this is not valid smarts %%")

    def test_empty_smiles_list_returns_empty_dataframe(self) -> None:
        """An empty input list should return an empty DataFrame."""
        df = decompose_rgroups([], BENZENE_CORE)

        assert isinstance(df, pd.DataFrame)
        assert df.empty

    def test_all_unmatched_returns_empty_dataframe(self) -> None:
        """When no compounds match the core, return an empty DataFrame."""
        df = decompose_rgroups([CYCLOHEXANE, "C1CCCC1"], BENZENE_CORE)

        assert isinstance(df, pd.DataFrame)
        assert df.empty

    def test_invalid_smiles_silently_excluded(self) -> None:
        """Invalid SMILES strings should be silently excluded."""
        smiles = [TOLUENE, "this_is_not_smiles", PHENOL]
        df = decompose_rgroups(smiles, BENZENE_CORE)

        # Only the two valid, matching molecules should appear
        assert len(df) == 2
        assert "this_is_not_smiles" not in df["smiles"].tolist()

    def test_core_column_present(self) -> None:
        """The Core column should be present in the output."""
        df = decompose_rgroups([TOLUENE, ANILINE], BENZENE_CORE)

        assert "Core" in df.columns

    def test_r1_values_are_strings(self) -> None:
        """R-group values should be SMILES strings (not Mol objects)."""
        df = decompose_rgroups([TOLUENE], BENZENE_CORE)

        assert not df.empty
        r1_val = df["R1"].iloc[0]
        assert isinstance(r1_val, str)


# ---------------------------------------------------------------------------
# rgroup_sar_table tests
# ---------------------------------------------------------------------------


class TestRgroupSarTable:
    """Tests for the rgroup_sar_table function."""

    def _make_compound(
        self,
        smiles: str,
        activity_value: float | None = None,
        activity_type: str = "IC50",
        units: str = "nM",
    ) -> Compound:
        """Helper: create a Compound with optional bioactivity."""
        c = Compound(smiles=smiles)
        if activity_value is not None:
            bio = Bioactivity(
                target_name="Test Target",
                activity_type=activity_type,
                value=activity_value,
                units=units,
            )
            c.bioactivities.append(bio)
        return c

    def test_sar_table_includes_activity_columns(self) -> None:
        """SAR table should contain activity_value, activity_units, pic50 columns."""
        compounds = [
            self._make_compound(TOLUENE, activity_value=100.0),
            self._make_compound(ANILINE, activity_value=500.0),
        ]
        df = rgroup_sar_table(compounds, BENZENE_CORE, activity_type="IC50")

        assert not df.empty
        assert "activity_value" in df.columns
        assert "activity_units" in df.columns
        assert "pic50" in df.columns

    def test_activity_values_joined_correctly(self) -> None:
        """Activity values should be correctly joined to the R-group data."""
        compounds = [
            self._make_compound(TOLUENE, activity_value=100.0),
            self._make_compound(PHENOL, activity_value=200.0),
        ]
        df = rgroup_sar_table(compounds, BENZENE_CORE, activity_type="IC50")

        toluene_row = df[df["smiles"] == TOLUENE]
        assert not toluene_row.empty
        assert toluene_row["activity_value"].iloc[0] == pytest.approx(100.0)

        phenol_row = df[df["smiles"] == PHENOL]
        assert not phenol_row.empty
        assert phenol_row["activity_value"].iloc[0] == pytest.approx(200.0)

    def test_pic50_computed_for_ic50(self) -> None:
        """pIC50 should be computed for IC50 entries with nM units."""
        import math

        compounds = [self._make_compound(TOLUENE, activity_value=100.0, units="nM")]
        df = rgroup_sar_table(compounds, BENZENE_CORE, activity_type="IC50")

        row = df[df["smiles"] == TOLUENE]
        assert not row.empty
        expected_pic50 = -math.log10(100.0 * 1e-9)
        assert row["pic50"].iloc[0] == pytest.approx(expected_pic50, abs=1e-6)

    def test_missing_activity_gives_nan(self) -> None:
        """Compounds without matching activity should have NaN activity values."""
        compounds = [
            self._make_compound(TOLUENE, activity_value=None),  # no activity
            self._make_compound(PHENOL, activity_value=50.0),
        ]
        df = rgroup_sar_table(compounds, BENZENE_CORE, activity_type="IC50")

        toluene_row = df[df["smiles"] == TOLUENE]
        assert not toluene_row.empty
        assert pd.isna(toluene_row["activity_value"].iloc[0])

    def test_activity_type_filter_respected(self) -> None:
        """Only activities matching the requested activity_type should be joined."""
        toluene = Compound(smiles=TOLUENE)
        bio_ki = Bioactivity(
            target_name="Target",
            activity_type="Ki",
            value=10.0,
            units="nM",
        )
        bio_ic50 = Bioactivity(
            target_name="Target",
            activity_type="IC50",
            value=999.0,
            units="nM",
        )
        toluene.bioactivities.extend([bio_ki, bio_ic50])

        df = rgroup_sar_table([toluene], BENZENE_CORE, activity_type="Ki")

        row = df[df["smiles"] == TOLUENE]
        assert not row.empty
        assert row["activity_value"].iloc[0] == pytest.approx(10.0)

    def test_empty_compounds_list_returns_empty_dataframe(self) -> None:
        """An empty compound list should return an empty DataFrame."""
        df = rgroup_sar_table([], BENZENE_CORE)

        assert isinstance(df, pd.DataFrame)
        assert df.empty

    def test_unmatched_compound_excluded_from_sar(self) -> None:
        """Compounds that do not match the core scaffold are excluded."""
        cyclohexane_compound = self._make_compound(CYCLOHEXANE, activity_value=1.0)
        toluene_compound = self._make_compound(TOLUENE, activity_value=100.0)

        df = rgroup_sar_table(
            [cyclohexane_compound, toluene_compound],
            BENZENE_CORE,
        )

        assert CYCLOHEXANE not in df["smiles"].tolist()
        assert TOLUENE in df["smiles"].tolist()


# ---------------------------------------------------------------------------
# CompoundCollection integration tests
# ---------------------------------------------------------------------------


class TestCompoundCollectionRgroup:
    """Tests for CompoundCollection.decompose_rgroups and sar_table methods."""

    def _make_collection(self, with_activity: bool = False) -> CompoundCollection:
        """Helper: build a test collection of substituted benzenes."""
        compounds = []
        for smi, val in [(TOLUENE, 100.0), (ANILINE, 200.0), (PHENOL, 50.0)]:
            c = Compound(smiles=smi)
            if with_activity:
                bio = Bioactivity(
                    target_name="Target",
                    activity_type="IC50",
                    value=val,
                    units="nM",
                )
                c.bioactivities.append(bio)
            compounds.append(c)
        return CompoundCollection(compounds=compounds)

    def test_collection_decompose_rgroups(self) -> None:
        """CompoundCollection.decompose_rgroups should delegate correctly."""
        collection = self._make_collection()
        df = collection.decompose_rgroups(BENZENE_CORE)

        assert not df.empty
        assert "R1" in df.columns
        assert len(df) == 3

    def test_collection_sar_table(self) -> None:
        """CompoundCollection.sar_table should include R-groups and activity."""
        collection = self._make_collection(with_activity=True)
        df = collection.sar_table(BENZENE_CORE, activity_type="IC50")

        assert not df.empty
        assert "R1" in df.columns
        assert "activity_value" in df.columns
        assert "pic50" in df.columns
        assert len(df) == 3

    def test_collection_sar_table_invalid_smarts_raises(self) -> None:
        """sar_table should propagate ValueError for invalid SMARTS."""
        collection = self._make_collection()
        with pytest.raises(ValueError, match="Invalid core SMARTS"):
            collection.sar_table("%%%invalid%%%")

    def test_collection_decompose_empty_collection(self) -> None:
        """An empty collection should return an empty DataFrame."""
        collection = CompoundCollection(compounds=[])
        df = collection.decompose_rgroups(BENZENE_CORE)

        assert isinstance(df, pd.DataFrame)
        assert df.empty
