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
    """Return a mock RDKit Mol object with realistic descriptor values.

    GetSubstructMatches returns [] by default so SMARTS-based checks
    (hERG basic N, AMES aromatic amine/nitroso) evaluate to False.
    """
    mol = MagicMock()
    mol.GetNumHeavyAtoms.return_value = 13
    mol.GetAtoms.return_value = [
        MagicMock(**{"GetIsAromatic.return_value": True, "GetAtomicNum.return_value": 6})
        for _ in range(6)
    ] + [
        MagicMock(**{"GetIsAromatic.return_value": False, "GetAtomicNum.return_value": 8})
        for _ in range(4)
    ]
    # Return empty tuple so all SMARTS substructure matches evaluate False
    mol.GetSubstructMatches.return_value = ()
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
            mock_chem.MolFromSmarts.return_value = MagicMock()
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            mock_desc.MolLogP.return_value = 1.2
            mock_desc.MolWt.return_value = 180.16
            mock_desc.TPSA.return_value = 63.6
            mock_rdmd.CalcNumRotatableBonds.return_value = 3
            mock_rdmd.CalcNumHBD.return_value = 1
            mock_rdmd.CalcNumHBA.return_value = 2
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
            mock_chem.MolFromSmarts.return_value = MagicMock()
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            mock_desc.MolLogP.return_value = 1.2
            mock_desc.MolWt.return_value = 180.16
            mock_desc.TPSA.return_value = 63.6
            mock_rdmd.CalcNumRotatableBonds.return_value = 3
            mock_rdmd.CalcNumHBD.return_value = 1
            mock_rdmd.CalcNumHBA.return_value = 2
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
            mock_chem.MolFromSmarts.return_value = MagicMock()
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            mock_desc.MolLogP.return_value = 1.2
            mock_desc.MolWt.return_value = 180.16
            mock_desc.TPSA.return_value = 63.6
            mock_rdmd.CalcNumRotatableBonds.return_value = 3
            mock_rdmd.CalcNumHBD.return_value = 1
            mock_rdmd.CalcNumHBA.return_value = 2
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
            mock_chem.MolFromSmarts.return_value = MagicMock()
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            mock_desc.MolLogP.return_value = 1.2
            mock_desc.MolWt.return_value = 180.16
            mock_desc.TPSA.return_value = 63.6
            mock_rdmd.CalcNumRotatableBonds.return_value = 3
            mock_rdmd.CalcNumHBD.return_value = 1
            mock_rdmd.CalcNumHBA.return_value = 2
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
            mock_chem.MolFromSmarts.return_value = MagicMock()
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            mock_admet_mod.ADMETModel.return_value = mock_model
            mock_desc.MolLogP.return_value = 1.2
            mock_desc.MolWt.return_value = 180.16
            mock_desc.TPSA.return_value = 63.6
            mock_rdmd.CalcNumRotatableBonds.return_value = 3
            mock_rdmd.CalcNumHBD.return_value = 1
            mock_rdmd.CalcNumHBA.return_value = 2
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
            mock_chem.MolFromSmarts.return_value = MagicMock()
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            mock_desc.MolLogP.return_value = 1.2
            mock_desc.MolWt.return_value = 180.16
            mock_desc.TPSA.return_value = 63.6
            mock_rdmd.CalcNumRotatableBonds.return_value = 3
            mock_rdmd.CalcNumHBD.return_value = 1
            mock_rdmd.CalcNumHBA.return_value = 2
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
            mock_rdmd.CalcNumHBA.return_value = 2
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
            mock_chem.MolFromSmarts.return_value = MagicMock()
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            mock_desc.MolLogP.return_value = 1.2
            mock_desc.MolWt.return_value = 180.16
            mock_desc.TPSA.return_value = 63.6
            mock_rdmd.CalcNumRotatableBonds.return_value = 3
            mock_rdmd.CalcNumHBD.return_value = 1
            mock_rdmd.CalcNumHBA.return_value = 2
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


# ---------------------------------------------------------------------------
# Tests: Scientific correctness (real RDKit molecules)
# Require RDKit to be installed; skipped otherwise.
# ---------------------------------------------------------------------------

try:
    from rdkit import Chem as _RDKitChem  # noqa: F401
    _RDKIT_AVAILABLE_FOR_TESTS = True
except ImportError:
    _RDKIT_AVAILABLE_FOR_TESTS = False

_requires_rdkit = pytest.mark.skipif(
    not _RDKIT_AVAILABLE_FOR_TESTS,
    reason="RDKit not installed; skipping scientific correctness tests",
)


@_requires_rdkit
class TestCYPIsoformDifferentiation:
    """Verify that CYP isoform functions score chemically distinct molecules differently.

    CYP1A2 risk drivers: planar aromatic scaffold (>= 3 rings), small-moderate MW.
    CYP2D6 risk drivers: basic aliphatic nitrogen, moderate MW and LogP.
    A molecule dominated by planarity should score higher on CYP1A2 than CYP2D6,
    and vice versa for a molecule with a basic amine but no extended planarity.
    """

    def test_cyp1a2_vs_cyp2d6_planar_aromatic(self):
        """Acridine (3 fused rings, no basic N) scores higher for CYP1A2 than CYP2D6."""
        from chemfuse.compute.admet import (
            _rule_cyp1a2_inhibition,
            _rule_cyp2d6_inhibition,
        )
        # Acridine: three linearly fused aromatic rings, no basic amine
        mol = _RDKitChem.MolFromSmiles("c1ccc2nc3ccccc3cc2c1")
        assert mol is not None, "Acridine SMILES failed to parse"

        cyp1a2 = _rule_cyp1a2_inhibition(mol)
        cyp2d6 = _rule_cyp2d6_inhibition(mol)

        # CYP1A2 should be elevated (planarity) and CYP2D6 low (no basic N)
        assert cyp1a2.value > cyp2d6.value, (
            f"Expected CYP1A2 risk ({cyp1a2.value}) > CYP2D6 risk ({cyp2d6.value}) "
            "for planar aromatic acridine"
        )

    def test_cyp2d6_vs_cyp1a2_basic_amine(self):
        """Diphenhydramine-like scaffold (basic N, 2 rings, LogP>3) favours CYP2D6 over CYP1A2.

        Molecule: CN(C)CCCC(c1ccccc1)c1ccccc1
        - 2 aromatic rings (< 3, so CYP1A2 planarity criterion fails)
        - LogP ~3.8 (> 3, so CYP1A2 low-LogP criterion also fails)
        - Basic tertiary N + MW in [200,500] + LogP in [1,4]: all CYP2D6 criteria met
        """
        from chemfuse.compute.admet import (
            _rule_cyp1a2_inhibition,
            _rule_cyp2d6_inhibition,
        )
        # Dimethylaminobutyl diphenylmethane: basic tertiary N, MW ~283, LogP ~3.8, 2 aryl rings
        mol = _RDKitChem.MolFromSmiles("CN(C)CCCC(c1ccccc1)c1ccccc1")
        assert mol is not None

        cyp2d6 = _rule_cyp2d6_inhibition(mol)
        cyp1a2 = _rule_cyp1a2_inhibition(mol)

        # CYP2D6 should be elevated (basic N + MW/LogP in range)
        # CYP1A2 should be lower (only 2 aromatic rings, LogP > 3)
        assert cyp2d6.value > cyp1a2.value, (
            f"Expected CYP2D6 risk ({cyp2d6.value}) > CYP1A2 risk ({cyp1a2.value}) "
            "for basic-amine-containing molecule with 2 aromatic rings and LogP > 3"
        )

    def test_cyp_confidence_range(self):
        """All CYP isoform predictions report confidence between 0.5 and 0.6."""
        from chemfuse.compute.admet import (
            _rule_cyp1a2_inhibition,
            _rule_cyp2c9_inhibition,
            _rule_cyp2c19_inhibition,
            _rule_cyp2d6_inhibition,
            _rule_cyp3a4_inhibition,
        )
        mol = _RDKitChem.MolFromSmiles(ASPIRIN_SMILES)
        assert mol is not None

        for fn in (
            _rule_cyp1a2_inhibition,
            _rule_cyp2c9_inhibition,
            _rule_cyp2c19_inhibition,
            _rule_cyp2d6_inhibition,
            _rule_cyp3a4_inhibition,
        ):
            pred = fn(mol)
            assert pred.confidence is not None
            assert 0.5 <= pred.confidence <= 0.6, (
                f"{fn.__name__} confidence {pred.confidence} outside [0.5, 0.6]"
            )


@_requires_rdkit
class TestAmesAromaticAmine:
    """Verify corrected aromatic amine detection in AMES mutagenicity rule.

    Ring nitrogen atoms (pyridine, pyrimidine) must NOT be flagged as aromatic
    amines.  An exocyclic -NH2/-NHR attached to an aromatic ring MUST be flagged.
    """

    def test_pyridine_does_not_trigger_aromatic_amine(self):
        """Pyridine ring N is not an aromatic amine; AMES score must not include amine contribution."""
        from chemfuse.compute.admet import _rule_ames

        mol = _RDKitChem.MolFromSmiles("c1ccncc1")  # pyridine
        assert mol is not None
        pred = _rule_ames(mol)
        # Pyridine has no nitro, no exocyclic amine, no nitroso -> score must be 0
        assert pred.value == 0.0, (
            f"Pyridine should have AMES score 0.0, got {pred.value}"
        )

    def test_aniline_triggers_aromatic_amine(self):
        """Aniline (-NH2 on benzene ring) must be flagged as an aromatic amine."""
        from chemfuse.compute.admet import _rule_ames

        mol = _RDKitChem.MolFromSmiles("c1ccc(N)cc1")  # aniline
        assert mol is not None
        pred = _rule_ames(mol)
        # Aromatic amine contributes 0.4 to score
        assert pred.value >= 0.4, (
            f"Aniline should have AMES score >= 0.4 (aromatic amine), got {pred.value}"
        )

    def test_pyrimidine_does_not_trigger_aromatic_amine(self):
        """Pyrimidine ring nitrogens must not be flagged as aromatic amines."""
        from chemfuse.compute.admet import _rule_ames

        mol = _RDKitChem.MolFromSmiles("c1cnccn1")  # pyrimidine
        assert mol is not None
        pred = _rule_ames(mol)
        assert pred.value == 0.0, (
            f"Pyrimidine should have AMES score 0.0, got {pred.value}"
        )

    def test_nitro_group_contributes_to_ames(self):
        """Nitrobenzene (nitro group) must score >= 0.5 on AMES."""
        from chemfuse.compute.admet import _rule_ames

        mol = _RDKitChem.MolFromSmiles("c1ccc([N+](=O)[O-])cc1")  # nitrobenzene
        assert mol is not None
        pred = _rule_ames(mol)
        assert pred.value >= 0.5, (
            f"Nitrobenzene should have AMES score >= 0.5 (nitro group), got {pred.value}"
        )

    def test_4_aminobiphenyl_triggers_aromatic_amine(self):
        """4-Aminobiphenyl (known mutagen) must be flagged via aromatic amine SMARTS."""
        from chemfuse.compute.admet import _rule_ames

        mol = _RDKitChem.MolFromSmiles("c1ccc(-c2ccc(N)cc2)cc1")  # 4-aminobiphenyl
        assert mol is not None
        pred = _rule_ames(mol)
        assert pred.value >= 0.4, (
            f"4-Aminobiphenyl should have AMES score >= 0.4, got {pred.value}"
        )


@_requires_rdkit
class TestHERGBasicNitrogen:
    """Verify that hERG liability uses basic nitrogen, not any nitrogen.

    Amide nitrogen must NOT trigger the basic-N risk factor.
    Piperidine (basic tertiary amine) MUST trigger the basic-N risk factor.
    """

    def test_amide_nitrogen_not_basic(self):
        """Amide N (CC(=O)NC) must not contribute to hERG basic-N risk."""
        from chemfuse.compute.admet import _rule_herg

        mol = _RDKitChem.MolFromSmiles("CC(=O)NC")  # N-methylacetamide
        assert mol is not None
        pred = _rule_herg(mol)
        # N-methylacetamide: low MW (~73), low LogP, no basic N -> risk should be 0
        assert pred.value == 0.0, (
            f"N-methylacetamide hERG risk should be 0.0 (amide N excluded), got {pred.value}"
        )

    def test_piperidine_triggers_basic_n(self):
        """Piperidine (basic secondary amine in ring) must contribute to hERG risk."""
        from chemfuse.compute.admet import _rule_herg

        mol = _RDKitChem.MolFromSmiles("C1CCNCC1")  # piperidine
        assert mol is not None
        pred = _rule_herg(mol)
        # Piperidine has basic N -> risk > 0 (basic_n flag is True)
        assert pred.value > 0.0, (
            f"Piperidine hERG risk should be > 0.0 (basic ring amine), got {pred.value}"
        )

    def test_sulfonamide_nitrogen_not_basic(self):
        """Sulfonamide N must not contribute to hERG basic-N risk."""
        from chemfuse.compute.admet import _rule_herg

        mol = _RDKitChem.MolFromSmiles("CS(=O)(=O)NC")  # N-methylmethanesulfonamide
        assert mol is not None
        pred = _rule_herg(mol)
        # Sulfonamide: low MW, no basic N -> risk should be 0
        assert pred.value == 0.0, (
            f"Sulfonamide N-methylmethanesulfonamide hERG risk should be 0.0, got {pred.value}"
        )

    def test_tertiary_amine_triggers_basic_n(self):
        """Trimethylamine (tertiary aliphatic amine) must trigger hERG basic-N risk."""
        from chemfuse.compute.admet import _rule_herg

        mol = _RDKitChem.MolFromSmiles("CN(C)C")  # trimethylamine
        assert mol is not None
        pred = _rule_herg(mol)
        # Trimethylamine has a basic tertiary N -> basic_n flag is True
        assert pred.value > 0.0, (
            f"Trimethylamine hERG risk should be > 0.0 (basic tertiary N), got {pred.value}"
        )

    def test_pyridine_nitrogen_not_basic_for_herg(self):
        """Pyridine aromatic ring N must not be counted as basic N for hERG."""
        from chemfuse.compute.admet import _rule_herg

        mol = _RDKitChem.MolFromSmiles("c1ccncc1")  # pyridine
        assert mol is not None
        pred = _rule_herg(mol)
        # Pyridine: low MW, moderate LogP, aromatic N excluded -> risk should be 0
        assert pred.value == 0.0, (
            f"Pyridine hERG risk should be 0.0 (aromatic ring N excluded), got {pred.value}"
        )
