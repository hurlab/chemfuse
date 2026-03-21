"""Tests for ADMET prediction module (SPEC-CF-005)."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from chemfuse.core.exceptions import ChemFuseError
from chemfuse.models.prediction import ADMETPrediction, ADMETProfile

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
CAFFEINE_SMILES = "Cn1cnc2c1c(=O)n(c(=O)n2C)C"
INVALID_SMILES = "THIS_IS_NOT_VALID"


# ---------------------------------------------------------------------------
# Helper: import admet module with mocked rdkit
# ---------------------------------------------------------------------------

def _make_mol_mock():
    """Return a mock RDKit Mol object with realistic descriptor values."""
    mol = MagicMock()
    mol.GetNumHeavyAtoms.return_value = 13
    mol.GetAtoms.return_value = [
        MagicMock(**{"GetIsAromatic.return_value": True, "GetAtomicNum.return_value": 6})
        for _ in range(6)
    ] + [
        MagicMock(**{"GetIsAromatic.return_value": False, "GetAtomicNum.return_value": 8})
        for _ in range(4)
    ]
    return mol


# ---------------------------------------------------------------------------
# Tests: rule-based predictions (no admet-ai, no RDKit mocked away)
# ---------------------------------------------------------------------------

class TestRuleBasedPredictions:
    """Tests for rule-based ADMET fallback (admet-ai not available)."""

    @patch("chemfuse.compute.admet._ADMET_AI_AVAILABLE", False)
    @patch("chemfuse.compute.admet._RDKIT_AVAILABLE", True)
    def test_predict_admet_returns_profile(self):
        """predict_admet returns an ADMETProfile for a valid SMILES."""
        from chemfuse.compute.admet import predict_admet

        with patch("chemfuse.compute.admet.Chem") as mock_chem, \
             patch("chemfuse.compute.admet.Descriptors") as mock_desc, \
             patch("chemfuse.compute.admet.rdMolDescriptors") as mock_rdmd:

            mol = _make_mol_mock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            mock_desc.MolLogP.return_value = 1.2
            mock_desc.MolWt.return_value = 180.16
            mock_desc.TPSA.return_value = 63.6
            mock_rdmd.CalcNumRotatableBonds.return_value = 3
            mock_rdmd.CalcNumHBD.return_value = 1
            mock_rdmd.CalcNumAromaticRings.return_value = 1

            result = predict_admet(ASPIRIN_SMILES)

        assert isinstance(result, ADMETProfile)
        assert result.smiles == ASPIRIN_SMILES
        assert len(result.predictions) > 0
        assert result.overall_score is not None

    @patch("chemfuse.compute.admet._ADMET_AI_AVAILABLE", False)
    @patch("chemfuse.compute.admet._RDKIT_AVAILABLE", True)
    def test_all_properties_present(self):
        """All 15 expected ADMET properties are populated."""
        from chemfuse.compute.admet import predict_admet

        expected_props = {
            "caco2_permeability", "hia", "bbb_permeability",
            "cyp1a2_inhibition", "cyp2c9_inhibition", "cyp2c19_inhibition",
            "cyp2d6_inhibition", "cyp3a4_inhibition", "half_life", "clearance",
            "herg_liability", "ames_mutagenicity", "dili",
            "solubility", "lipophilicity",
        }

        with patch("chemfuse.compute.admet.Chem") as mock_chem, \
             patch("chemfuse.compute.admet.Descriptors") as mock_desc, \
             patch("chemfuse.compute.admet.rdMolDescriptors") as mock_rdmd:

            mol = _make_mol_mock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            mock_desc.MolLogP.return_value = 1.2
            mock_desc.MolWt.return_value = 180.16
            mock_desc.TPSA.return_value = 63.6
            mock_rdmd.CalcNumRotatableBonds.return_value = 3
            mock_rdmd.CalcNumHBD.return_value = 1
            mock_rdmd.CalcNumAromaticRings.return_value = 1

            result = predict_admet(ASPIRIN_SMILES)

        assert expected_props.issubset(set(result.predictions.keys()))

    @patch("chemfuse.compute.admet._ADMET_AI_AVAILABLE", False)
    @patch("chemfuse.compute.admet._RDKIT_AVAILABLE", True)
    def test_all_predictions_are_rule_based(self):
        """All predictions have method='rule-based' when admet-ai is absent."""
        from chemfuse.compute.admet import predict_admet

        with patch("chemfuse.compute.admet.Chem") as mock_chem, \
             patch("chemfuse.compute.admet.Descriptors") as mock_desc, \
             patch("chemfuse.compute.admet.rdMolDescriptors") as mock_rdmd:

            mol = _make_mol_mock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            mock_desc.MolLogP.return_value = 1.2
            mock_desc.MolWt.return_value = 180.16
            mock_desc.TPSA.return_value = 63.6
            mock_rdmd.CalcNumRotatableBonds.return_value = 3
            mock_rdmd.CalcNumHBD.return_value = 1
            mock_rdmd.CalcNumAromaticRings.return_value = 1

            result = predict_admet(ASPIRIN_SMILES)

        for pred in result.predictions.values():
            assert pred.method == "rule-based"

    @patch("chemfuse.compute.admet._ADMET_AI_AVAILABLE", False)
    @patch("chemfuse.compute.admet._RDKIT_AVAILABLE", True)
    def test_invalid_smiles_raises(self):
        """Invalid SMILES raises ChemFuseError."""
        from chemfuse.compute.admet import predict_admet

        with patch("chemfuse.compute.admet.Chem") as mock_chem:
            mock_chem.MolFromSmiles.return_value = None
            with pytest.raises(ChemFuseError, match="Invalid SMILES"):
                predict_admet(INVALID_SMILES)

    @patch("chemfuse.compute.admet._ADMET_AI_AVAILABLE", False)
    @patch("chemfuse.compute.admet._RDKIT_AVAILABLE", True)
    def test_empty_smiles_raises(self):
        """Empty SMILES string raises ChemFuseError."""
        from chemfuse.compute.admet import predict_admet

        with pytest.raises(ChemFuseError, match="empty"):
            predict_admet("")

    @patch("chemfuse.compute.admet._ADMET_AI_AVAILABLE", False)
    @patch("chemfuse.compute.admet._RDKIT_AVAILABLE", True)
    def test_prediction_has_confidence(self):
        """Each prediction includes a confidence score."""
        from chemfuse.compute.admet import predict_admet

        with patch("chemfuse.compute.admet.Chem") as mock_chem, \
             patch("chemfuse.compute.admet.Descriptors") as mock_desc, \
             patch("chemfuse.compute.admet.rdMolDescriptors") as mock_rdmd:

            mol = _make_mol_mock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            mock_desc.MolLogP.return_value = 1.2
            mock_desc.MolWt.return_value = 180.16
            mock_desc.TPSA.return_value = 63.6
            mock_rdmd.CalcNumRotatableBonds.return_value = 3
            mock_rdmd.CalcNumHBD.return_value = 1
            mock_rdmd.CalcNumAromaticRings.return_value = 1

            result = predict_admet(ASPIRIN_SMILES)

        for pred in result.predictions.values():
            assert pred.confidence is not None
            assert 0.0 <= pred.confidence <= 1.0


# ---------------------------------------------------------------------------
# Tests: ML-based predictions (admet-ai available)
# ---------------------------------------------------------------------------

class TestMLBasedPredictions:
    """Tests for ML ADMET path (admet-ai mocked as available)."""

    def test_all_predictions_are_ml(self):
        """When admet-ai is available, all predictions have method='ml'."""
        from chemfuse.compute.admet import predict_admet

        mock_model = MagicMock()
        mock_model.predict.return_value = {
            "caco2_permeability": 12.5,
            "bbb_permeability": 0.3,
            "solubility": -3.1,
            "hia": 85.0,
        }

        with patch("chemfuse.compute.admet._ADMET_AI_AVAILABLE", True), \
             patch("chemfuse.compute.admet._RDKIT_AVAILABLE", True), \
             patch("chemfuse.compute.admet.Chem") as mock_chem, \
             patch("chemfuse.compute.admet.admet_ai") as mock_admet_mod:

            mock_chem.MolFromSmiles.return_value = _make_mol_mock()
            mock_admet_mod.ADMETModel.return_value = mock_model

            result = predict_admet(ASPIRIN_SMILES)

        for pred in result.predictions.values():
            assert pred.method == "ml"

    def test_ml_failure_falls_back_to_rule_based(self):
        """If admet-ai model.predict fails, falls back to rule-based."""
        from chemfuse.compute.admet import predict_admet

        mock_model = MagicMock()
        mock_model.predict.side_effect = RuntimeError("model exploded")

        with patch("chemfuse.compute.admet._ADMET_AI_AVAILABLE", True), \
             patch("chemfuse.compute.admet._RDKIT_AVAILABLE", True), \
             patch("chemfuse.compute.admet.Chem") as mock_chem, \
             patch("chemfuse.compute.admet.Descriptors") as mock_desc, \
             patch("chemfuse.compute.admet.rdMolDescriptors") as mock_rdmd, \
             patch("chemfuse.compute.admet.admet_ai") as mock_admet_mod:

            mol = _make_mol_mock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            mock_admet_mod.ADMETModel.return_value = mock_model
            mock_desc.MolLogP.return_value = 1.2
            mock_desc.MolWt.return_value = 180.16
            mock_desc.TPSA.return_value = 63.6
            mock_rdmd.CalcNumRotatableBonds.return_value = 3
            mock_rdmd.CalcNumHBD.return_value = 1
            mock_rdmd.CalcNumAromaticRings.return_value = 1

            result = predict_admet(ASPIRIN_SMILES)

        # Should still return a profile (rule-based fallback)
        assert isinstance(result, ADMETProfile)
        for pred in result.predictions.values():
            assert pred.method == "rule-based"


# ---------------------------------------------------------------------------
# Tests: Batch predictions
# ---------------------------------------------------------------------------

class TestBatchPredictions:
    """Tests for predict_admet_batch."""

    @patch("chemfuse.compute.admet._ADMET_AI_AVAILABLE", False)
    @patch("chemfuse.compute.admet._RDKIT_AVAILABLE", True)
    def test_batch_returns_correct_count(self):
        """predict_admet_batch returns one profile per input SMILES."""
        from chemfuse.compute.admet import predict_admet_batch

        smiles_list = [ASPIRIN_SMILES, CAFFEINE_SMILES, ASPIRIN_SMILES]

        with patch("chemfuse.compute.admet.Chem") as mock_chem, \
             patch("chemfuse.compute.admet.Descriptors") as mock_desc, \
             patch("chemfuse.compute.admet.rdMolDescriptors") as mock_rdmd:

            mol = _make_mol_mock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            mock_desc.MolLogP.return_value = 1.2
            mock_desc.MolWt.return_value = 180.16
            mock_desc.TPSA.return_value = 63.6
            mock_rdmd.CalcNumRotatableBonds.return_value = 3
            mock_rdmd.CalcNumHBD.return_value = 1
            mock_rdmd.CalcNumAromaticRings.return_value = 1

            results = predict_admet_batch(smiles_list)

        assert len(results) == 3
        for r in results:
            assert isinstance(r, ADMETProfile)

    @patch("chemfuse.compute.admet._ADMET_AI_AVAILABLE", False)
    @patch("chemfuse.compute.admet._RDKIT_AVAILABLE", True)
    def test_batch_handles_invalid_smiles_gracefully(self):
        """predict_admet_batch does not raise on invalid SMILES; returns empty profile."""
        from chemfuse.compute.admet import predict_admet_batch

        def side_effect(smi):
            return None if smi == INVALID_SMILES else _make_mol_mock()

        with patch("chemfuse.compute.admet.Chem") as mock_chem, \
             patch("chemfuse.compute.admet.Descriptors") as mock_desc, \
             patch("chemfuse.compute.admet.rdMolDescriptors") as mock_rdmd:

            mock_chem.MolFromSmiles.side_effect = side_effect
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            mock_desc.MolLogP.return_value = 1.2
            mock_desc.MolWt.return_value = 180.16
            mock_desc.TPSA.return_value = 63.6
            mock_rdmd.CalcNumRotatableBonds.return_value = 3
            mock_rdmd.CalcNumHBD.return_value = 1
            mock_rdmd.CalcNumAromaticRings.return_value = 1

            results = predict_admet_batch([ASPIRIN_SMILES, INVALID_SMILES])

        assert len(results) == 2
        # Second result is the fallback empty profile
        assert results[1].smiles == INVALID_SMILES
        assert results[1].predictions == {}

    @patch("chemfuse.compute.admet._ADMET_AI_AVAILABLE", False)
    @patch("chemfuse.compute.admet._RDKIT_AVAILABLE", True)
    def test_batch_100_smiles(self):
        """predict_admet_batch handles 100 SMILES input."""
        from chemfuse.compute.admet import predict_admet_batch

        with patch("chemfuse.compute.admet.Chem") as mock_chem, \
             patch("chemfuse.compute.admet.Descriptors") as mock_desc, \
             patch("chemfuse.compute.admet.rdMolDescriptors") as mock_rdmd:

            mol = _make_mol_mock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            mock_desc.MolLogP.return_value = 1.2
            mock_desc.MolWt.return_value = 180.16
            mock_desc.TPSA.return_value = 63.6
            mock_rdmd.CalcNumRotatableBonds.return_value = 3
            mock_rdmd.CalcNumHBD.return_value = 1
            mock_rdmd.CalcNumAromaticRings.return_value = 1

            results = predict_admet_batch([ASPIRIN_SMILES] * 100)

        assert len(results) == 100


# ---------------------------------------------------------------------------
# Tests: ADMETPrediction and ADMETProfile models
# ---------------------------------------------------------------------------

class TestADMETModels:
    """Tests for Pydantic model definitions."""

    def test_admet_prediction_model(self):
        """ADMETPrediction model accepts all expected fields."""
        pred = ADMETPrediction(
            property_name="bbb_permeability",
            value=0.75,
            unit=None,
            confidence=0.6,
            method="rule-based",
            category="medium",
        )
        assert pred.property_name == "bbb_permeability"
        assert pred.method == "rule-based"
        assert pred.category == "medium"

    def test_admet_prediction_default_method(self):
        """ADMETPrediction defaults method to 'rule-based'."""
        pred = ADMETPrediction(property_name="solubility", value=-3.0)
        assert pred.method == "rule-based"

    def test_admet_profile_model(self):
        """ADMETProfile model serializes correctly."""
        profile = ADMETProfile(
            smiles=ASPIRIN_SMILES,
            predictions={
                "solubility": ADMETPrediction(
                    property_name="solubility", value=-2.5, method="ml"
                )
            },
            overall_score=0.8,
        )
        assert profile.smiles == ASPIRIN_SMILES
        assert "solubility" in profile.predictions
        assert profile.overall_score == 0.8

    def test_admet_profile_empty_predictions(self):
        """ADMETProfile allows empty predictions dict."""
        profile = ADMETProfile(smiles=ASPIRIN_SMILES)
        assert profile.predictions == {}
        assert profile.overall_score is None

    def test_admet_profile_ml_method_label(self):
        """ADMETPrediction can be marked as 'ml'."""
        pred = ADMETPrediction(
            property_name="hia", value=90.0, method="ml", confidence=0.95
        )
        assert pred.method == "ml"
