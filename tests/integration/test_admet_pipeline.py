"""Integration tests for ADMET prediction pipeline (SPEC-CF-005)."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from chemfuse.models.collection import CompoundCollection
from chemfuse.models.compound import Compound, CompoundProperties
from chemfuse.models.prediction import ADMETPrediction, ADMETProfile

ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
CAFFEINE_SMILES = "Cn1cnc2c1c(=O)n(c(=O)n2C)C"
IBUPROFEN_SMILES = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"


def _make_compound(smiles: str, name: str) -> Compound:
    return Compound(
        smiles=smiles,
        name=name,
        sources=["test"],
        properties=CompoundProperties(molecular_weight=180.0),
    )


def _make_collection(*pairs: tuple[str, str]) -> CompoundCollection:
    compounds = [_make_compound(s, n) for s, n in pairs]
    return CompoundCollection(compounds=compounds)


def _make_mol_mock() -> MagicMock:
    """Return a MagicMock RDKit mol with numeric descriptor side-effects."""
    mol = MagicMock()
    mol.GetNumHeavyAtoms.return_value = 13
    mol.GetAtoms.return_value = []
    mol.HasSubstructMatch.return_value = False
    mol.GetSubstructMatches.return_value = ()
    return mol


def _mock_descriptors(mock_desc: MagicMock, mock_rdmd: MagicMock) -> None:
    """Configure descriptor mocks to return realistic float values."""
    mock_desc.MolLogP.return_value = 1.2
    mock_desc.MolWt.return_value = 180.16
    mock_desc.ExactMolWt.return_value = 180.16
    mock_desc.TPSA.return_value = 63.6
    mock_desc.NumHDonors.return_value = 1
    mock_desc.NumHAcceptors.return_value = 4
    mock_rdmd.CalcNumRotatableBonds.return_value = 3
    mock_rdmd.CalcNumAromaticRings.return_value = 1
    mock_rdmd.CalcNumHBD.return_value = 1
    mock_rdmd.CalcNumHBA.return_value = 4


# ---------------------------------------------------------------------------
# ADMETPrediction / ADMETProfile model tests
# ---------------------------------------------------------------------------

class TestADMETPredictionModel:
    def test_valid_prediction_creation(self) -> None:
        """ADMETPrediction model can be instantiated with required fields."""
        pred = ADMETPrediction(property_name="hia", value=0.85, method="rule-based")
        assert pred.property_name == "hia"
        assert pred.value == pytest.approx(0.85)
        assert pred.method == "rule-based"

    def test_prediction_with_all_fields(self) -> None:
        """ADMETPrediction stores all optional fields."""
        pred = ADMETPrediction(
            property_name="caco2_permeability",
            value=12.5,
            unit="nm/s",
            confidence=0.9,
            method="ml",
            category="absorption",
        )
        assert pred.unit == "nm/s"
        assert pred.confidence == pytest.approx(0.9)
        assert pred.category == "absorption"


class TestADMETProfileModel:
    def test_empty_profile_creation(self) -> None:
        """ADMETProfile can be created with just smiles."""
        profile = ADMETProfile(smiles=ASPIRIN_SMILES)
        assert profile.smiles == ASPIRIN_SMILES
        assert profile.predictions == {}
        assert profile.overall_score is None

    def test_profile_with_predictions(self) -> None:
        """ADMETProfile stores dict of ADMETPrediction objects."""
        pred = ADMETPrediction(property_name="hia", value=0.9, method="rule-based")
        profile = ADMETProfile(smiles=ASPIRIN_SMILES, predictions={"hia": pred}, overall_score=0.75)
        assert "hia" in profile.predictions
        assert profile.overall_score == pytest.approx(0.75)


# ---------------------------------------------------------------------------
# predict_admet function tests
# ---------------------------------------------------------------------------

class TestPredictAdmet:
    def test_returns_admet_profile(self) -> None:
        """predict_admet returns an ADMETProfile for a valid SMILES."""
        with patch("chemfuse.compute.admet._RDKIT_AVAILABLE", True), \
             patch("chemfuse.compute.admet.Chem") as mock_chem, \
             patch("chemfuse.compute.admet.Descriptors") as mock_desc, \
             patch("chemfuse.compute.admet.rdMolDescriptors") as mock_rdmd:

            mol = _make_mol_mock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            _mock_descriptors(mock_desc, mock_rdmd)

            from chemfuse.compute.admet import predict_admet
            profile = predict_admet(ASPIRIN_SMILES)

        assert isinstance(profile, ADMETProfile)
        assert profile.smiles == ASPIRIN_SMILES

    def test_profile_contains_predictions(self) -> None:
        """predict_admet result has predictions dict with expected keys."""
        with patch("chemfuse.compute.admet._RDKIT_AVAILABLE", True), \
             patch("chemfuse.compute.admet.Chem") as mock_chem, \
             patch("chemfuse.compute.admet.Descriptors") as mock_desc, \
             patch("chemfuse.compute.admet.rdMolDescriptors") as mock_rdmd:

            mol = _make_mol_mock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            _mock_descriptors(mock_desc, mock_rdmd)

            from chemfuse.compute.admet import predict_admet
            profile = predict_admet(ASPIRIN_SMILES)

        assert len(profile.predictions) > 0
        expected_keys = {"hia", "bbb_permeability", "solubility"}
        assert any(k in profile.predictions for k in expected_keys)

    def test_overall_score_is_float(self) -> None:
        """predict_admet overall_score is a float in [0, 1]."""
        with patch("chemfuse.compute.admet._RDKIT_AVAILABLE", True), \
             patch("chemfuse.compute.admet.Chem") as mock_chem, \
             patch("chemfuse.compute.admet.Descriptors") as mock_desc, \
             patch("chemfuse.compute.admet.rdMolDescriptors") as mock_rdmd:

            mol = _make_mol_mock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            _mock_descriptors(mock_desc, mock_rdmd)

            from chemfuse.compute.admet import predict_admet
            profile = predict_admet(ASPIRIN_SMILES)

        if profile.overall_score is not None:
            assert 0.0 <= profile.overall_score <= 1.0

    def test_invalid_smiles_raises(self) -> None:
        """predict_admet raises ChemFuseError for invalid SMILES."""
        with patch("chemfuse.compute.admet._RDKIT_AVAILABLE", True), \
             patch("chemfuse.compute.admet.Chem") as mock_chem:

            mock_chem.MolFromSmiles.return_value = None

            from chemfuse.compute.admet import predict_admet
            from chemfuse.core.exceptions import ChemFuseError

            with pytest.raises(ChemFuseError):
                predict_admet("BADSMILES!!!")

    def test_empty_smiles_raises(self) -> None:
        """predict_admet raises ChemFuseError for empty SMILES."""
        from chemfuse.compute.admet import predict_admet
        from chemfuse.core.exceptions import ChemFuseError

        with pytest.raises(ChemFuseError):
            predict_admet("")

    def test_without_rdkit_returns_placeholder_profile(self) -> None:
        """predict_admet without RDKit returns a profile with empty predictions."""
        with patch("chemfuse.compute.admet._RDKIT_AVAILABLE", False), \
             patch("chemfuse.compute.admet._ADMET_AI_AVAILABLE", False):

            from chemfuse.compute.admet import predict_admet
            profile = predict_admet(ASPIRIN_SMILES)

        assert isinstance(profile, ADMETProfile)
        assert profile.smiles == ASPIRIN_SMILES


class TestPredictAdmetBatch:
    def test_returns_list_of_profiles(self) -> None:
        """predict_admet_batch returns one profile per input SMILES."""
        smiles_list = [ASPIRIN_SMILES, CAFFEINE_SMILES, IBUPROFEN_SMILES]

        with patch("chemfuse.compute.admet._RDKIT_AVAILABLE", True), \
             patch("chemfuse.compute.admet.Chem") as mock_chem, \
             patch("chemfuse.compute.admet.Descriptors") as mock_desc, \
             patch("chemfuse.compute.admet.rdMolDescriptors") as mock_rdmd:

            mol = _make_mol_mock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            _mock_descriptors(mock_desc, mock_rdmd)

            from chemfuse.compute.admet import predict_admet_batch
            profiles = predict_admet_batch(smiles_list)

        assert len(profiles) == 3
        assert all(isinstance(p, ADMETProfile) for p in profiles)

    def test_invalid_smiles_in_batch_returns_empty_profile(self) -> None:
        """predict_admet_batch handles invalid SMILES gracefully."""
        with patch("chemfuse.compute.admet._RDKIT_AVAILABLE", True), \
             patch("chemfuse.compute.admet.Chem") as mock_chem, \
             patch("chemfuse.compute.admet.Descriptors") as mock_desc, \
             patch("chemfuse.compute.admet.rdMolDescriptors") as mock_rdmd:

            mol = _make_mol_mock()
            # First call (aspirin) returns valid mol; second call (invalid) returns None
            mock_chem.MolFromSmiles.side_effect = [mol, None]
            mock_chem.MolToSmiles.return_value = ASPIRIN_SMILES
            _mock_descriptors(mock_desc, mock_rdmd)

            from chemfuse.compute.admet import predict_admet_batch
            profiles = predict_admet_batch([ASPIRIN_SMILES, "INVALIDSMILES"])

        assert len(profiles) == 2
        assert profiles[1].smiles == "INVALIDSMILES"

    def test_empty_list_returns_empty(self) -> None:
        """predict_admet_batch with empty input returns empty list."""
        from chemfuse.compute.admet import predict_admet_batch
        assert predict_admet_batch([]) == []


# ---------------------------------------------------------------------------
# CompoundCollection.predict_admet() integration
# ---------------------------------------------------------------------------

class TestCollectionPredictAdmet:
    def test_predict_admet_stores_on_compounds(self) -> None:
        """CompoundCollection.predict_admet stores profiles on each compound."""
        collection = _make_collection(
            (ASPIRIN_SMILES, "aspirin"),
            (CAFFEINE_SMILES, "caffeine"),
        )

        mock_profile = ADMETProfile(
            smiles=ASPIRIN_SMILES,
            predictions={"hia": ADMETPrediction(property_name="hia", value=0.9, method="rule-based")},
            overall_score=0.75,
        )

        with patch("chemfuse.compute.admet.predict_admet", return_value=mock_profile):
            collection.predict_admet()

        # Verify profile is stored on compound instances
        for compound in collection:
            assert hasattr(compound, "_admet_profile")
            profile = compound._admet_profile
            assert isinstance(profile, ADMETProfile)

    def test_predict_admet_skips_empty_smiles(self) -> None:
        """predict_admet skips compounds without SMILES."""
        from chemfuse.models.compound import Compound

        compound_no_smiles = Compound(smiles="", name="empty", sources=["test"])
        collection = CompoundCollection(compounds=[compound_no_smiles])

        with patch("chemfuse.compute.admet.predict_admet") as mock_predict:
            collection.predict_admet()

        mock_predict.assert_not_called()

    def test_predict_admet_handles_errors_gracefully(self) -> None:
        """predict_admet continues despite individual compound failures."""
        collection = _make_collection(
            (ASPIRIN_SMILES, "aspirin"),
            (CAFFEINE_SMILES, "caffeine"),
        )

        from chemfuse.core.exceptions import ChemFuseError

        def _raise_for_caffeine(smiles: str) -> ADMETProfile:
            if "caffeine" in smiles or smiles == CAFFEINE_SMILES:
                raise ChemFuseError("Prediction failed")
            return ADMETProfile(smiles=smiles)

        with patch("chemfuse.compute.admet.predict_admet", side_effect=_raise_for_caffeine):
            # Should not raise
            collection.predict_admet()

        # aspirin should have profile
        aspirin = collection.compounds[0]
        assert hasattr(aspirin, "_admet_profile")

        # caffeine should NOT have profile (error was swallowed)
        caffeine = collection.compounds[1]
        assert not hasattr(caffeine, "_admet_profile")
