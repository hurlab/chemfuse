"""Tests for Bioactivity and BindingMeasurement models."""

from __future__ import annotations

import pytest

from chemfuse.models.bioactivity import (
    RECOGNIZED_ACTIVITY_TYPES,
    BindingMeasurement,
    Bioactivity,
    normalize_to_nm,
)
from chemfuse.models.compound import Compound

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


class TestNormalizeToNm:
    """Tests for the normalize_to_nm helper function."""

    def test_nm_stays_same(self) -> None:
        """100 nM stays 100."""
        assert normalize_to_nm(100.0, "nM") == pytest.approx(100.0)

    def test_um_to_nm(self) -> None:
        """0.1 uM becomes 100 nM."""
        assert normalize_to_nm(0.1, "uM") == pytest.approx(100.0)

    def test_m_to_nm(self) -> None:
        """1e-7 M becomes 100 nM."""
        assert normalize_to_nm(1e-7, "M") == pytest.approx(100.0)

    def test_mm_to_nm(self) -> None:
        """0.0001 mM becomes 100 nM."""
        assert normalize_to_nm(0.0001, "mM") == pytest.approx(100.0)

    def test_pm_to_nm(self) -> None:
        """100000 pM becomes 100 nM."""
        assert normalize_to_nm(100000.0, "pM") == pytest.approx(100.0)

    def test_micromolar_string(self) -> None:
        """'micromolar' is treated as uM."""
        assert normalize_to_nm(0.1, "micromolar") == pytest.approx(100.0)

    def test_nanomolar_string(self) -> None:
        """'nanomolar' is treated as nM."""
        assert normalize_to_nm(100.0, "nanomolar") == pytest.approx(100.0)

    def test_none_units_returns_none(self) -> None:
        """None units returns None."""
        assert normalize_to_nm(100.0, None) is None

    def test_none_value_returns_none(self) -> None:
        """None value returns None."""
        assert normalize_to_nm(None, "nM") is None

    def test_unrecognized_units_returns_none(self) -> None:
        """Unrecognized units returns None."""
        assert normalize_to_nm(100.0, "ug/mL") is None

    def test_case_insensitive_units(self) -> None:
        """Unit matching is case-insensitive."""
        assert normalize_to_nm(0.1, "UM") == pytest.approx(100.0)
        assert normalize_to_nm(0.1, "UМ") is None  # non-ASCII M stays None

    @pytest.mark.parametrize("value,units,expected", [
        (100.0, "nM", 100.0),
        (0.1, "uM", 100.0),
        (1e-7, "M", 100.0),
        (0.0001, "mM", 100.0),
        (100000.0, "pM", 100.0),
    ])
    def test_normalize_parametrized(
        self,
        value: float,
        units: str,
        expected: float,
    ) -> None:
        """Parametrized unit conversion all produce 100 nM."""
        assert normalize_to_nm(value, units) == pytest.approx(expected)


class TestBioactivityNormalization:
    """Tests for auto-normalization in Bioactivity model."""

    def test_ic50_100nm_pic50_is_7(self) -> None:
        """IC50 of 100 nM -> pIC50 == 7.0."""
        ba = Bioactivity(
            target_name="Target",
            activity_type="IC50",
            value=100.0,
            units="nM",
        )
        assert ba.value_nm == pytest.approx(100.0)
        assert ba.pic50 == pytest.approx(7.0)
        assert ba.pki is None
        assert ba.is_normalized is True

    def test_ki_10nm_pki_is_8(self) -> None:
        """Ki of 10 nM -> pKi == 8.0."""
        ba = Bioactivity(
            target_name="Target",
            activity_type="Ki",
            value=10.0,
            units="nM",
        )
        assert ba.value_nm == pytest.approx(10.0)
        assert ba.pki == pytest.approx(8.0)
        assert ba.pic50 is None
        assert ba.is_normalized is True

    def test_ec50_does_not_set_pic50_or_pki(self) -> None:
        """EC50 with nM units: value_nm is set but pic50 and pki are None."""
        ba = Bioactivity(
            target_name="Target",
            activity_type="EC50",
            value=100.0,
            units="nM",
        )
        assert ba.value_nm == pytest.approx(100.0)
        assert ba.pic50 is None
        assert ba.pki is None

    def test_no_units_no_normalization(self) -> None:
        """No units means no normalization."""
        ba = Bioactivity(
            target_name="Target",
            activity_type="IC50",
            value=100.0,
            units=None,
        )
        assert ba.value_nm is None
        assert ba.pic50 is None
        assert ba.is_normalized is False

    def test_uM_ic50_normalization(self) -> None:
        """IC50 of 0.1 uM -> value_nm == 100 nM, pIC50 == 7.0."""
        ba = Bioactivity(
            target_name="Target",
            activity_type="IC50",
            value=0.1,
            units="uM",
        )
        assert ba.value_nm == pytest.approx(100.0)
        assert ba.pic50 == pytest.approx(7.0)

    def test_unrecognized_units_no_normalization(self) -> None:
        """Unrecognized units leaves value_nm None."""
        ba = Bioactivity(
            target_name="Target",
            activity_type="IC50",
            value=1.0,
            units="ug/mL",
        )
        assert ba.value_nm is None
        assert ba.is_normalized is False

    def test_new_fields_in_model_dump(self) -> None:
        """New fields appear in model_dump output."""
        ba = Bioactivity(target_name="T", activity_type="Ki")
        d = ba.model_dump()
        assert "value_nm" in d
        assert "pic50" in d
        assert "pki" in d
        assert "is_normalized" in d


class TestBestActivity:
    """Tests for Compound.best_activity() method."""

    def _make_compound_with_activities(self) -> Compound:
        """Create a compound with mixed bioactivity data."""
        from chemfuse.models.bioactivity import Bioactivity

        c = Compound(smiles="CC(=O)Oc1ccccc1C(O)=O", name="Aspirin")
        c.bioactivities = [
            Bioactivity(
                target_name="COX-1",
                activity_type="IC50",
                value=1670.0,
                units="nM",
            ),
            Bioactivity(
                target_name="COX-2",
                activity_type="IC50",
                value=278000.0,
                units="nM",
            ),
            Bioactivity(
                target_name="COX-1",
                activity_type="Ki",
                value=500.0,
                units="nM",
            ),
        ]
        return c

    def test_best_activity_returns_lowest_nm(self) -> None:
        """best_activity returns the lowest value_nm for the activity type."""
        c = self._make_compound_with_activities()
        result = c.best_activity(activity_type="IC50")
        assert result == pytest.approx(1670.0)

    def test_best_activity_filters_by_target(self) -> None:
        """best_activity filters by target name substring."""
        c = self._make_compound_with_activities()
        result = c.best_activity(target="COX-2", activity_type="IC50")
        assert result == pytest.approx(278000.0)

    def test_best_activity_ki_type(self) -> None:
        """best_activity works for Ki type."""
        c = self._make_compound_with_activities()
        result = c.best_activity(activity_type="Ki")
        assert result == pytest.approx(500.0)

    def test_best_activity_no_match_returns_none(self) -> None:
        """best_activity returns None when no matching activities."""
        c = self._make_compound_with_activities()
        result = c.best_activity(target="p38", activity_type="IC50")
        assert result is None

    def test_best_activity_empty_bioactivities(self) -> None:
        """best_activity returns None for compound with no bioactivities."""
        c = Compound(smiles="C", name="methane")
        assert c.best_activity() is None

    def test_best_activity_no_normalized_values(self) -> None:
        """best_activity returns None when no value_nm available."""
        from chemfuse.models.bioactivity import Bioactivity

        c = Compound(smiles="C", name="methane")
        c.bioactivities = [
            Bioactivity(
                target_name="Target",
                activity_type="IC50",
                value=100.0,
                units=None,  # No units -> no normalization
            ),
        ]
        assert c.best_activity() is None


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
