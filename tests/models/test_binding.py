"""Additional tests for BindingMeasurement model edge cases."""

from __future__ import annotations

import pytest

from chemfuse.models.bioactivity import BindingMeasurement


class TestBindingMeasurementEdgeCases:
    """Edge case tests for BindingMeasurement."""

    def test_uniprot_stored_correctly(self) -> None:
        """target_uniprot is stored and retrieved correctly."""
        bm = BindingMeasurement(target_name="COX-1", target_uniprot="P23219")
        assert bm.target_uniprot == "P23219"

    def test_reference_stored_correctly(self) -> None:
        """reference field is stored correctly."""
        bm = BindingMeasurement(target_name="T", reference="10.1234/test")
        assert bm.reference == "10.1234/test"

    def test_all_relation_fields(self) -> None:
        """All relation fields can be set."""
        bm = BindingMeasurement(
            target_name="T",
            ki=1000.0,
            ki_relation=">",
            kd=100.0,
            kd_relation="<",
            ic50=50.0,
            ic50_relation=">=",
            ec50=200.0,
            ec50_relation="<=",
        )
        assert bm.ki_relation == ">"
        assert bm.kd_relation == "<"
        assert bm.ic50_relation == ">="
        assert bm.ec50_relation == "<="

    def test_best_affinity_excludes_high_inequality(self) -> None:
        """best_affinity uses numeric value regardless of relation flag."""
        # The numeric value is still used for comparison
        bm = BindingMeasurement(
            target_name="T",
            ki=10000.0,
            ki_relation=">",  # means ">10000", but we store numeric 10000
            ic50=100.0,
        )
        # best_affinity returns min of numeric values
        assert bm.best_affinity == pytest.approx(100.0)

    def test_source_can_be_overridden(self) -> None:
        """source field can be set to non-default values."""
        bm = BindingMeasurement(target_name="T", source="chembl")
        assert bm.source == "chembl"

    def test_units_fields_optional(self) -> None:
        """Unit fields are all optional."""
        bm = BindingMeasurement(target_name="T", ki=1.0)
        assert bm.ki_units is None
        assert bm.kd_units is None
        assert bm.ic50_units is None
        assert bm.ec50_units is None

    def test_units_fields_can_be_set(self) -> None:
        """Unit fields can be set."""
        bm = BindingMeasurement(
            target_name="T",
            ki=1.0,
            ki_units="nM",
            ic50=10.0,
            ic50_units="uM",
        )
        assert bm.ki_units == "nM"
        assert bm.ic50_units == "uM"
