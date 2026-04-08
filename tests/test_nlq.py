"""Tests for chemfuse.nlq — Natural Language Query interface."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from chemfuse.models.collection import CompoundCollection
from chemfuse.models.compound import Compound

# Common SMILES used across tests
ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
ETHANOL_SMILES = "CCO"
ACETIC_ACID_SMILES = "CC(=O)O"


def _make_collection(*smiles_list: str) -> CompoundCollection:
    """Build a minimal CompoundCollection from SMILES strings."""
    compounds = [
        Compound(smiles=s, name=s, sources=["pubchem"])
        for s in smiles_list
    ]
    return CompoundCollection(compounds=compounds, query="test", sources=["pubchem"])


class TestAskSearch:
    """Tests for 'search for X' / 'find X' patterns."""

    def test_search_for_aspirin_calls_search(self):
        """'search for aspirin' calls chemfuse.search() with name query."""
        expected = _make_collection(ASPIRIN_SMILES)
        with patch("chemfuse.search", return_value=expected) as mock_search:
            from chemfuse.nlq import ask
            result = ask("search for aspirin")

        mock_search.assert_called_once_with("aspirin", query_type="name")
        assert result is expected

    def test_find_compound_calls_search(self):
        """'find caffeine' also routes to search()."""
        expected = _make_collection()
        with patch("chemfuse.search", return_value=expected) as mock_search:
            from chemfuse.nlq import ask
            result = ask("find caffeine")

        mock_search.assert_called_once_with("caffeine", query_type="name")
        assert result is expected

    def test_search_with_smiles_uses_smiles_query_type(self):
        """'search for CC(=O)O' uses query_type='smiles'."""
        expected = _make_collection(ACETIC_ACID_SMILES)
        with patch("chemfuse.search", return_value=expected) as mock_search:
            from chemfuse.nlq import ask
            result = ask(f"search for {ACETIC_ACID_SMILES}")

        mock_search.assert_called_once_with(ACETIC_ACID_SMILES, query_type="smiles")
        assert result is expected


class TestAskFindSimilar:
    """Tests for 'similar to X' / 'find compounds like X' patterns."""

    def test_similar_to_smiles_calls_find_similar(self):
        """'find compounds similar to CC(=O)O' calls find_similar()."""
        expected = _make_collection(ACETIC_ACID_SMILES)
        with patch("chemfuse.find_similar", return_value=expected) as mock_sim:
            from chemfuse.nlq import ask
            result = ask(f"find compounds similar to {ACETIC_ACID_SMILES}")

        mock_sim.assert_called_once_with(ACETIC_ACID_SMILES)
        assert result is expected

    def test_compounds_like_smiles(self):
        """'compounds like CC(=O)O' pattern works."""
        expected = _make_collection()
        with patch("chemfuse.find_similar", return_value=expected) as mock_sim:
            from chemfuse.nlq import ask
            ask(f"compounds like {ACETIC_ACID_SMILES}")

        mock_sim.assert_called_once_with(ACETIC_ACID_SMILES)

    def test_similar_to_aspirin_smiles(self):
        """Full aspirin SMILES in a similarity query is routed correctly."""
        expected = _make_collection()
        with patch("chemfuse.find_similar", return_value=expected) as mock_sim:
            from chemfuse.nlq import ask
            ask(f"similar to {ASPIRIN_SMILES}")

        mock_sim.assert_called_once_with(ASPIRIN_SMILES)


class TestAskAdmet:
    """Tests for 'ADMET for X' / 'predict ADMET X' patterns."""

    def test_admet_for_smiles_calls_predict_admet(self):
        """'ADMET for CC(=O)Oc1ccccc1C(=O)O' calls predict_admet."""
        mock_profile = MagicMock()
        with patch("chemfuse.compute.admet.predict_admet", return_value=mock_profile) as mock_admet:
            from chemfuse.nlq import ask
            result = ask(f"ADMET for {ASPIRIN_SMILES}")

        mock_admet.assert_called_once_with(ASPIRIN_SMILES)
        assert result is mock_profile

    def test_predict_admet_prefix(self):
        """'predict ADMET CC(=O)O' is also accepted."""
        mock_profile = MagicMock()
        with patch("chemfuse.compute.admet.predict_admet", return_value=mock_profile) as mock_admet:
            from chemfuse.nlq import ask
            ask(f"predict ADMET {ACETIC_ACID_SMILES}")

        mock_admet.assert_called_once_with(ACETIC_ACID_SMILES)

    def test_admet_for_smiles_acetic_acid(self):
        """'ADMET for CC(=O)O' correctly extracts acetic acid SMILES."""
        mock_profile = MagicMock()
        with patch("chemfuse.compute.admet.predict_admet", return_value=mock_profile) as mock_admet:
            from chemfuse.nlq import ask
            ask(f"ADMET for {ACETIC_ACID_SMILES}")

        mock_admet.assert_called_once_with(ACETIC_ACID_SMILES)


class TestAskDrugLikeness:
    """Tests for 'is X drug-like?' / 'drug-likeness of X' patterns."""

    def test_is_smiles_drug_like_calls_check(self):
        """'is CC(=O)O drug-like' calls check_drug_likeness."""
        mock_dl = MagicMock()
        with patch("chemfuse.compute.druglikeness.check_drug_likeness", return_value=mock_dl), \
             patch("chemfuse.compute.descriptors.compute_descriptors", return_value={}):
            from chemfuse.nlq import ask
            result = ask(f"is {ACETIC_ACID_SMILES} drug-like")

        assert result is mock_dl

    def test_drug_likeness_of_smiles(self):
        """'drug-likeness of CC(=O)O' pattern works."""
        mock_dl = MagicMock()
        with patch("chemfuse.compute.druglikeness.check_drug_likeness", return_value=mock_dl), \
             patch("chemfuse.compute.descriptors.compute_descriptors", return_value={}):
            from chemfuse.nlq import ask
            result = ask(f"drug-likeness of {ACETIC_ACID_SMILES}")

        assert result is mock_dl

    def test_check_druglikeness_pattern(self):
        """'check drug-likeness for CC(=O)O' works."""
        mock_dl = MagicMock()
        with patch("chemfuse.compute.druglikeness.check_drug_likeness", return_value=mock_dl), \
             patch("chemfuse.compute.descriptors.compute_descriptors", return_value={}):
            from chemfuse.nlq import ask
            result = ask(f"check drug-likeness for {ACETIC_ACID_SMILES}")

        assert result is mock_dl


class TestAskInvalidQuery:
    """Tests for queries that cannot be parsed."""

    def test_empty_string_raises(self):
        """Empty query raises ValueError."""
        from chemfuse.nlq import ask
        with pytest.raises(ValueError):
            ask("")

    def test_whitespace_only_raises(self):
        """Whitespace-only query raises ValueError."""
        from chemfuse.nlq import ask
        with pytest.raises(ValueError):
            ask("   ")

    def test_unrecognized_query_raises_value_error(self):
        """A query that matches no pattern raises ValueError."""
        from chemfuse.nlq import ask
        with pytest.raises(ValueError, match="Cannot parse query"):
            ask("zzzzz xyzzy frobulate")

    def test_error_message_contains_examples(self):
        """ValueError message contains usage hints."""
        from chemfuse.nlq import ask
        with pytest.raises(ValueError) as exc_info:
            ask("zzzzz xyzzy frobulate")
        assert "search" in str(exc_info.value).lower() or "admet" in str(exc_info.value).lower()


class TestAskCompare:
    """Tests for 'compare X and Y' pattern."""

    def test_compare_two_smiles_returns_dict(self):
        """'compare CC(=O)O and CCO' returns a dict with both compounds."""
        with patch("chemfuse.compute.descriptors.compute_descriptors", return_value={"MolWt": 60.0}), \
             patch("chemfuse.compute.druglikeness.check_drug_likeness", return_value=MagicMock()):
            from chemfuse.nlq import ask
            result = ask(f"compare {ACETIC_ACID_SMILES} and {ETHANOL_SMILES}")

        assert isinstance(result, dict)
        assert result["compound_a"] == ACETIC_ACID_SMILES
        assert result["compound_b"] == ETHANOL_SMILES

    def test_compare_result_has_properties_keys(self):
        """Comparison dict includes properties_a and properties_b."""
        with patch("chemfuse.compute.descriptors.compute_descriptors", return_value={}), \
             patch("chemfuse.compute.druglikeness.check_drug_likeness", return_value=MagicMock()):
            from chemfuse.nlq import ask
            result = ask(f"compare {ACETIC_ACID_SMILES} and {ETHANOL_SMILES}")

        assert "properties_a" in result
        assert "properties_b" in result

    def test_compare_result_has_druglikeness_keys(self):
        """Comparison dict includes druglikeness_a and druglikeness_b."""
        mock_dl = MagicMock()
        with patch("chemfuse.compute.descriptors.compute_descriptors", return_value={}), \
             patch("chemfuse.compute.druglikeness.check_drug_likeness", return_value=mock_dl):
            from chemfuse.nlq import ask
            result = ask(f"compare {ACETIC_ACID_SMILES} and {ETHANOL_SMILES}")

        assert "druglikeness_a" in result
        assert "druglikeness_b" in result


class TestAskPublicApi:
    """Tests verifying 'ask' is exported in the public API."""

    def test_ask_importable_from_chemfuse(self):
        """ask is importable directly from chemfuse."""
        import chemfuse
        assert hasattr(chemfuse, "ask")
        assert callable(chemfuse.ask)

    def test_ask_in_all(self):
        """ask is in chemfuse.__all__."""
        import chemfuse
        assert "ask" in chemfuse.__all__


class TestIsSmiles:
    """Unit tests for the internal SMILES detection heuristic."""

    def test_smiles_with_parentheses_detected(self):
        """CC(=O)O is detected as SMILES."""
        from chemfuse.nlq import _is_smiles
        assert _is_smiles("CC(=O)O") is True

    def test_smiles_with_brackets_detected(self):
        """[NH4+] is detected as SMILES."""
        from chemfuse.nlq import _is_smiles
        assert _is_smiles("[NH4+]") is True

    def test_plain_name_not_smiles(self):
        """'aspirin' is not detected as SMILES."""
        from chemfuse.nlq import _is_smiles
        assert _is_smiles("aspirin") is False

    def test_empty_string_not_smiles(self):
        """Empty string is not SMILES."""
        from chemfuse.nlq import _is_smiles
        assert _is_smiles("") is False

    def test_aromatic_smiles_detected(self):
        """c1ccccc1 (benzene) is detected as SMILES via lowercase atoms."""
        from chemfuse.nlq import _is_smiles
        assert _is_smiles("c1ccccc1") is True
