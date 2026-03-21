"""Tests for Bioactivity and BindingMeasurement models."""

from __future__ import annotations

import pytest

from chemfuse.models.bioactivity import (
    RECOGNIZED_ACTIVITY_TYPES,
    BindingMeasurement,
    Bioactivity,
)

# --- Bioactivity tests ---


class TestBioactivity:
    """Tests for the Bioactivity Pydantic model."""

    def test_basic_creation(self) -> None:
        """Create a Bioactivity with required fields."""
        ba = Bioactivity(
            target_name="Cyclooxygenase-1",
            activity_type="IC50",
            value=1.5,
            units="uM",
        )
        assert ba.target_name == "Cyclooxygenase-1"
        assert ba.activity_type == "IC50"
        assert ba.value == pytest.approx(1.5)
        assert ba.units == "uM"

    def test_optional_fields_default_none(self) -> None:
        """Optional fields default to None."""
        ba = Bioactivity(target_name="Target", activity_type="Ki")
        assert ba.target_id is None
        assert ba.value is None
        assert ba.units is None
        assert ba.relation is None
        assert ba.assay_type is None
        assert ba.reference is None

    def test_source_defaults_to_empty_string(self) -> None:
        """Source defaults to empty string."""
        ba = Bioactivity(target_name="T", activity_type="Ki")
        assert ba.source == ""

    def test_activity_type_normalization_recognized(self) -> None:
        """Recognized activity types are normalized to canonical form."""
        for atype in ["IC50", "Ki", "Kd", "EC50", "ED50", "GI50"]:
            ba = Bioactivity(target_name="T", activity_type=atype)
            assert ba.activity_type == atype

    def test_activity_type_case_insensitive(self) -> None:
        """Activity type normalization is case-insensitive."""
        ba = Bioactivity(target_name="T", activity_type="ic50")
        assert ba.activity_type == "IC50"

        ba2 = Bioactivity(target_name="T", activity_type="KI")
        assert ba2.activity_type == "Ki"

    def test_unknown_activity_type_becomes_other(self) -> None:
        """Unrecognized activity types are normalized to 'other'."""
        ba = Bioactivity(target_name="T", activity_type="MIC")
        assert ba.activity_type == "other"

    def test_other_activity_type_explicit(self) -> None:
        """Explicitly setting 'other' works."""
        ba = Bioactivity(target_name="T", activity_type="other")
        assert ba.activity_type == "other"

    def test_serialization_to_dict(self) -> None:
        """to_dict returns a flat dictionary."""
        ba = Bioactivity(
            target_name="COX-1",
            target_id="CHEMBL247",
            activity_type="IC50",
            value=10.0,
            units="nM",
            relation="=",
            assay_type="Binding",
            source="chembl",
            reference="CHEMBL123",
        )
        d = ba.to_dict()
        assert isinstance(d, dict)
        assert d["target_name"] == "COX-1"
        assert d["activity_type"] == "IC50"
        assert d["value"] == pytest.approx(10.0)
        assert d["units"] == "nM"
        assert d["source"] == "chembl"

    def test_model_dump_includes_all_fields(self) -> None:
        """model_dump includes all expected keys."""
        ba = Bioactivity(target_name="T", activity_type="Ki")
        d = ba.model_dump()
        expected_keys = {
            "target_name", "target_id", "activity_type", "value", "units",
            "relation", "assay_type", "source", "reference",
        }
        assert expected_keys.issubset(set(d.keys()))

    def test_all_recognized_types_in_constant(self) -> None:
        """RECOGNIZED_ACTIVITY_TYPES contains the expected values."""
        expected = {"IC50", "Ki", "Kd", "EC50", "ED50", "GI50", "other"}
        assert expected == RECOGNIZED_ACTIVITY_TYPES


class TestBindingMeasurement:
    """Tests for the BindingMeasurement Pydantic model."""

    def test_basic_creation(self) -> None:
        """Create a BindingMeasurement with required fields."""
        bm = BindingMeasurement(
            target_name="Cyclooxygenase-1",
            ki=1670.0,
        )
        assert bm.target_name == "Cyclooxygenase-1"
        assert bm.ki == pytest.approx(1670.0)

    def test_all_affinity_fields_optional(self) -> None:
        """All affinity fields are optional."""
        bm = BindingMeasurement(target_name="T")
        assert bm.ki is None
        assert bm.kd is None
        assert bm.ic50 is None
        assert bm.ec50 is None

    def test_source_defaults_to_bindingdb(self) -> None:
        """Source defaults to 'bindingdb'."""
        bm = BindingMeasurement(target_name="T")
        assert bm.source == "bindingdb"

    def test_best_affinity_returns_minimum(self) -> None:
        """best_affinity returns the smallest non-None value."""
        bm = BindingMeasurement(
            target_name="T",
            ki=1000.0,
            kd=500.0,
            ic50=250.0,
            ec50=750.0,
        )
        assert bm.best_affinity == pytest.approx(250.0)

    def test_best_affinity_ki_only(self) -> None:
        """best_affinity works with only Ki."""
        bm = BindingMeasurement(target_name="T", ki=1670.0)
        assert bm.best_affinity == pytest.approx(1670.0)

    def test_best_affinity_kd_wins_when_lower(self) -> None:
        """best_affinity returns Kd when it is lower than Ki."""
        bm = BindingMeasurement(target_name="T", ki=500.0, kd=50.0)
        assert bm.best_affinity == pytest.approx(50.0)

    def test_best_affinity_none_when_no_values(self) -> None:
        """best_affinity returns None when no measurements are set."""
        bm = BindingMeasurement(target_name="T")
        assert bm.best_affinity is None

    def test_best_affinity_ignores_none(self) -> None:
        """best_affinity ignores None values."""
        bm = BindingMeasurement(target_name="T", ki=None, kd=100.0, ic50=None)
        assert bm.best_affinity == pytest.approx(100.0)

    def test_serialization_to_dict(self) -> None:
        """to_dict returns a flat dictionary."""
        bm = BindingMeasurement(
            target_name="COX-1",
            target_uniprot="P23219",
            ki=1670.0,
            ic50=278000.0,
            source="bindingdb",
        )
        d = bm.to_dict()
        assert isinstance(d, dict)
        assert d["target_name"] == "COX-1"
        assert d["target_uniprot"] == "P23219"
        assert d["ki"] == pytest.approx(1670.0)
        assert d["ic50"] == pytest.approx(278000.0)
        assert d["source"] == "bindingdb"

    def test_model_dump_includes_all_fields(self) -> None:
        """model_dump includes all expected keys."""
        bm = BindingMeasurement(target_name="T")
        d = bm.model_dump()
        expected_keys = {
            "target_name", "target_uniprot", "ki", "kd", "ic50", "ec50",
            "ki_units", "kd_units", "ic50_units", "ec50_units",
            "ki_relation", "kd_relation", "ic50_relation", "ec50_relation",
            "source", "reference",
        }
        assert expected_keys.issubset(set(d.keys()))

    def test_inequality_fields(self) -> None:
        """Inequality relation fields are stored correctly."""
        bm = BindingMeasurement(
            target_name="T",
            ki=10000.0,
            ki_relation=">",
        )
        assert bm.ki_relation == ">"

    @pytest.mark.parametrize("ki,kd,ic50,ec50,expected", [
        (10.0, 50.0, 100.0, 200.0, 10.0),
        (None, 50.0, 100.0, 200.0, 50.0),
        (None, None, 100.0, 200.0, 100.0),
        (None, None, None, 200.0, 200.0),
        (None, None, None, None, None),
        (1.0, 1.0, 1.0, 1.0, 1.0),  # all equal
    ])
    def test_best_affinity_parametrized(
        self,
        ki: float | None,
        kd: float | None,
        ic50: float | None,
        ec50: float | None,
        expected: float | None,
    ) -> None:
        """best_affinity returns the minimum of all non-None values."""
        bm = BindingMeasurement(target_name="T", ki=ki, kd=kd, ic50=ic50, ec50=ec50)
        if expected is None:
            assert bm.best_affinity is None
        else:
            assert bm.best_affinity == pytest.approx(expected)
