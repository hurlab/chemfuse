"""Tests for chemfuse.compute.synthetic (SA score, Fsp3, NP-likeness)."""

from __future__ import annotations

from unittest.mock import patch

import pytest

from chemfuse.compute.synthetic import fraction_sp3, np_likeness, synthetic_accessibility
from chemfuse.core.exceptions import OptionalDependencyError

# ---------------------------------------------------------------------------
# Molecule SMILES used across tests
# ---------------------------------------------------------------------------

ETHANOL = "CCO"
BENZENE = "c1ccccc1"
CYCLOHEXANE = "C1CCCCC1"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"

# Taxol (paclitaxel) - complex, difficult to synthesize
TAXOL = (
    "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(=C2OC1=O)C)(C)C)OC(=O)C5=CC=CC=C5)(CO4)OC(=O)C)"
    "C)OC(=O)[C@@H]([C@@H](C6=CC=CC=C6)NC(=O)C7=CC=CC=C7)O"
)

INVALID_SMILES = "not_a_molecule"
EMPTY_SMILES = ""


# ---------------------------------------------------------------------------
# synthetic_accessibility tests
# ---------------------------------------------------------------------------


class TestSyntheticAccessibility:
    def test_ethanol_is_easy(self) -> None:
        """Ethanol (CCO) should have SA score well below 3.0."""
        score = synthetic_accessibility(ETHANOL)
        assert score is not None
        assert score < 3.0, f"Expected SA < 3.0 for ethanol, got {score}"

    def test_taxol_is_hard(self) -> None:
        """Taxol-like complex molecule should have SA score > 5.0."""
        score = synthetic_accessibility(TAXOL)
        assert score is not None
        assert score > 5.0, f"Expected SA > 5.0 for Taxol, got {score}"

    def test_score_range(self) -> None:
        """SA score must be in [1.0, 10.0]."""
        score = synthetic_accessibility(ASPIRIN)
        assert score is not None
        assert 1.0 <= score <= 10.0

    def test_invalid_smiles_returns_none(self) -> None:
        """Invalid SMILES should return None, not raise."""
        result = synthetic_accessibility(INVALID_SMILES)
        assert result is None

    def test_empty_smiles_returns_none(self) -> None:
        """Empty SMILES string should return None."""
        result = synthetic_accessibility(EMPTY_SMILES)
        assert result is None

    def test_raises_optional_dependency_error_without_rdkit(self) -> None:
        """OptionalDependencyError raised when RDKit is not available."""
        with patch("chemfuse.compute.synthetic._RDKIT_AVAILABLE", False):
            with pytest.raises(OptionalDependencyError):
                synthetic_accessibility(ETHANOL)


# ---------------------------------------------------------------------------
# fraction_sp3 tests
# ---------------------------------------------------------------------------


class TestFractionSp3:
    def test_cyclohexane_fsp3_is_one(self) -> None:
        """Cyclohexane has only sp3 carbons -> Fsp3 = 1.0."""
        result = fraction_sp3(CYCLOHEXANE)
        assert result is not None
        assert result == pytest.approx(1.0), f"Expected Fsp3=1.0 for cyclohexane, got {result}"

    def test_benzene_fsp3_is_zero(self) -> None:
        """Benzene has only aromatic carbons -> Fsp3 = 0.0."""
        result = fraction_sp3(BENZENE)
        assert result is not None
        assert result == pytest.approx(0.0), f"Expected Fsp3=0.0 for benzene, got {result}"

    def test_aspirin_fsp3_is_low(self) -> None:
        """Aspirin is mostly aromatic -> Fsp3 should be low (< 0.2)."""
        result = fraction_sp3(ASPIRIN)
        assert result is not None
        # Aspirin has 1 sp3 carbon (acetyl methyl) out of 9 carbons total
        assert result < 0.2, f"Expected low Fsp3 for aspirin, got {result}"

    def test_ethanol_fsp3_is_one(self) -> None:
        """Ethanol has only sp3 carbons -> Fsp3 = 1.0."""
        result = fraction_sp3(ETHANOL)
        assert result is not None
        assert result == pytest.approx(1.0)

    def test_result_in_range(self) -> None:
        """Fsp3 must be in [0.0, 1.0]."""
        result = fraction_sp3(ASPIRIN)
        assert result is not None
        assert 0.0 <= result <= 1.0

    def test_invalid_smiles_returns_none(self) -> None:
        """Invalid SMILES should return None."""
        result = fraction_sp3(INVALID_SMILES)
        assert result is None

    def test_empty_smiles_returns_none(self) -> None:
        """Empty SMILES should return None."""
        result = fraction_sp3(EMPTY_SMILES)
        assert result is None

    def test_raises_optional_dependency_error_without_rdkit(self) -> None:
        """OptionalDependencyError raised when RDKit is not available."""
        with patch("chemfuse.compute.synthetic._RDKIT_AVAILABLE", False):
            with pytest.raises(OptionalDependencyError):
                fraction_sp3(ETHANOL)


# ---------------------------------------------------------------------------
# np_likeness tests
# ---------------------------------------------------------------------------


class TestNpLikeness:
    def test_returns_float_or_none(self) -> None:
        """NP-likeness returns a float or None (both are valid)."""
        result = np_likeness(ASPIRIN)
        assert result is None or isinstance(result, float)

    def test_invalid_smiles_returns_none(self) -> None:
        """Invalid SMILES should return None."""
        result = np_likeness(INVALID_SMILES)
        assert result is None

    def test_empty_smiles_returns_none(self) -> None:
        """Empty SMILES should return None."""
        result = np_likeness(EMPTY_SMILES)
        assert result is None

    def test_raises_optional_dependency_error_without_rdkit(self) -> None:
        """OptionalDependencyError raised when RDKit is not available."""
        with patch("chemfuse.compute.synthetic._RDKIT_AVAILABLE", False):
            with pytest.raises(OptionalDependencyError):
                np_likeness(ETHANOL)

    def test_valid_smiles_with_npscorer_unavailable_returns_none(self) -> None:
        """When NP_Score module is unavailable, returns None silently."""
        with (
            patch("chemfuse.compute.synthetic._npscorer", None),
            patch("chemfuse.compute.synthetic._np_model", None),
            patch("chemfuse.compute.synthetic._npscorer_loaded", True),
        ):
            result = np_likeness(ETHANOL)
            assert result is None


# ---------------------------------------------------------------------------
# CompoundProperties model field tests
# ---------------------------------------------------------------------------


class TestCompoundPropertiesFields:
    def test_sa_score_field_exists(self) -> None:
        """CompoundProperties has sa_score field defaulting to None."""
        from chemfuse.models.compound import CompoundProperties

        props = CompoundProperties()
        assert props.sa_score is None

    def test_fsp3_field_exists(self) -> None:
        """CompoundProperties has fsp3 field defaulting to None."""
        from chemfuse.models.compound import CompoundProperties

        props = CompoundProperties()
        assert props.fsp3 is None

    def test_sa_score_can_be_set(self) -> None:
        """CompoundProperties accepts a float sa_score value."""
        from chemfuse.models.compound import CompoundProperties

        props = CompoundProperties(sa_score=2.5)
        assert props.sa_score == pytest.approx(2.5)

    def test_fsp3_can_be_set(self) -> None:
        """CompoundProperties accepts a float fsp3 value."""
        from chemfuse.models.compound import CompoundProperties

        props = CompoundProperties(fsp3=0.75)
        assert props.fsp3 == pytest.approx(0.75)
