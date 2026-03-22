"""Tests for chemfuse.compute.descriptors."""

from __future__ import annotations

import sys
from unittest.mock import patch

import pytest

ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
CAFFEINE_SMILES = "Cn1c(=O)c2c(ncn2C)n(c1=O)C"
INVALID_SMILES = "not_a_smiles_string!!!"


class TestComputeDescriptors:
    """Tests for compute_descriptors() with RDKit available."""

    def test_returns_dict(self):
        from chemfuse.compute.descriptors import compute_descriptors
        result = compute_descriptors(ASPIRIN_SMILES)
        assert isinstance(result, dict)

    def test_returns_200_plus_descriptors(self):
        from chemfuse.compute.descriptors import compute_descriptors
        result = compute_descriptors(ASPIRIN_SMILES)
        assert len(result) >= 200, f"Expected 200+ descriptors, got {len(result)}"

    def test_values_are_floats(self):
        from chemfuse.compute.descriptors import compute_descriptors
        result = compute_descriptors(ASPIRIN_SMILES)
        for name, value in result.items():
            assert isinstance(value, (int, float)), f"Descriptor {name!r} is not numeric"

    def test_caffeine_returns_dict(self):
        from chemfuse.compute.descriptors import compute_descriptors
        result = compute_descriptors(CAFFEINE_SMILES)
        assert isinstance(result, dict)
        assert len(result) >= 200

    def test_invalid_smiles_returns_empty_dict(self):
        from chemfuse.compute.descriptors import compute_descriptors
        result = compute_descriptors(INVALID_SMILES)
        assert result == {}

    def test_invalid_smiles_logs_warning(self, caplog):
        import logging

        from chemfuse.compute.descriptors import compute_descriptors
        with caplog.at_level(logging.WARNING, logger="chemfuse.compute.descriptors"):
            compute_descriptors(INVALID_SMILES)
        assert any("Invalid SMILES" in r.message for r in caplog.records)

    def test_rdkit_not_installed_raises_import_error(self):
        with patch.dict(sys.modules, {"rdkit": None, "rdkit.Chem": None}):
            # Re-import with mocked module to simulate missing rdkit
            import chemfuse.compute.descriptors as desc_mod
            original = desc_mod.RDKIT_AVAILABLE
            desc_mod.RDKIT_AVAILABLE = False
            try:
                with pytest.raises(ImportError, match="pip install chemfuse"):
                    desc_mod.compute_descriptors(ASPIRIN_SMILES)
            finally:
                desc_mod.RDKIT_AVAILABLE = original


class TestComputeDescriptorsBatch:
    """Tests for compute_descriptors_batch()."""

    def test_batch_returns_list(self):
        from chemfuse.compute.descriptors import compute_descriptors_batch
        results = compute_descriptors_batch([ASPIRIN_SMILES, CAFFEINE_SMILES])
        assert isinstance(results, list)
        assert len(results) == 2

    def test_batch_each_element_is_dict(self):
        from chemfuse.compute.descriptors import compute_descriptors_batch
        results = compute_descriptors_batch([ASPIRIN_SMILES, CAFFEINE_SMILES])
        for r in results:
            assert isinstance(r, dict)

    def test_batch_invalid_smiles_yields_empty_dict(self):
        from chemfuse.compute.descriptors import compute_descriptors_batch
        results = compute_descriptors_batch([ASPIRIN_SMILES, INVALID_SMILES])
        # First should have descriptors, second should be empty
        assert len(results[0]) >= 200
        assert results[1] == {}

    def test_batch_single_element(self):
        from chemfuse.compute.descriptors import compute_descriptors_batch
        results = compute_descriptors_batch([ASPIRIN_SMILES])
        assert len(results) == 1
        assert len(results[0]) >= 200

    def test_batch_empty_list(self):
        from chemfuse.compute.descriptors import compute_descriptors_batch
        results = compute_descriptors_batch([])
        assert results == []

    def test_batch_rdkit_not_installed_raises(self):
        import chemfuse.compute.descriptors as desc_mod
        original = desc_mod.RDKIT_AVAILABLE
        desc_mod.RDKIT_AVAILABLE = False
        try:
            with pytest.raises(ImportError, match="pip install chemfuse"):
                desc_mod.compute_descriptors_batch([ASPIRIN_SMILES])
        finally:
            desc_mod.RDKIT_AVAILABLE = original


class TestIsValidSmiles:
    """Tests for is_valid_smiles()."""

    def test_valid_aspirin(self):
        from chemfuse.compute.descriptors import is_valid_smiles
        assert is_valid_smiles(ASPIRIN_SMILES) is True

    def test_valid_methane(self):
        from chemfuse.compute.descriptors import is_valid_smiles
        assert is_valid_smiles("C") is True

    def test_invalid_smiles_returns_false(self):
        from chemfuse.compute.descriptors import is_valid_smiles
        assert is_valid_smiles(INVALID_SMILES) is False

    def test_rdkit_not_available_returns_false(self):
        """When RDKit is not available, cannot validate — return False."""
        import chemfuse.compute.descriptors as desc_mod
        original = desc_mod.RDKIT_AVAILABLE
        desc_mod.RDKIT_AVAILABLE = False
        try:
            assert desc_mod.is_valid_smiles(INVALID_SMILES) is False
        finally:
            desc_mod.RDKIT_AVAILABLE = original


class TestSmilesToInchi:
    """Tests for smiles_to_inchi()."""

    def test_aspirin_returns_inchi(self):
        from chemfuse.compute.descriptors import smiles_to_inchi
        result = smiles_to_inchi(ASPIRIN_SMILES)
        assert result is not None
        assert result.startswith("InChI=")

    def test_invalid_smiles_returns_none(self):
        from chemfuse.compute.descriptors import smiles_to_inchi
        result = smiles_to_inchi(INVALID_SMILES)
        assert result is None

    def test_rdkit_not_installed_raises(self):
        import chemfuse.compute.descriptors as desc_mod
        original = desc_mod.RDKIT_AVAILABLE
        desc_mod.RDKIT_AVAILABLE = False
        try:
            with pytest.raises(ImportError, match="pip install chemfuse"):
                desc_mod.smiles_to_inchi(ASPIRIN_SMILES)
        finally:
            desc_mod.RDKIT_AVAILABLE = original


class TestSmilesToInchikey:
    """Tests for smiles_to_inchikey()."""

    def test_aspirin_returns_inchikey(self):
        from chemfuse.compute.descriptors import smiles_to_inchikey
        result = smiles_to_inchikey(ASPIRIN_SMILES)
        assert result is not None
        # InChIKey format: XXXXXXXXXXXXXX-XXXXXXXXXX-X
        assert len(result) == 27
        assert result.count("-") == 2

    def test_known_aspirin_inchikey(self):
        from chemfuse.compute.descriptors import smiles_to_inchikey
        result = smiles_to_inchikey(ASPIRIN_SMILES)
        assert result == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

    def test_invalid_smiles_returns_none(self):
        from chemfuse.compute.descriptors import smiles_to_inchikey
        result = smiles_to_inchikey(INVALID_SMILES)
        assert result is None

    def test_rdkit_not_installed_raises(self):
        import chemfuse.compute.descriptors as desc_mod
        original = desc_mod.RDKIT_AVAILABLE
        desc_mod.RDKIT_AVAILABLE = False
        try:
            with pytest.raises(ImportError, match="pip install chemfuse"):
                desc_mod.smiles_to_inchikey(ASPIRIN_SMILES)
        finally:
            desc_mod.RDKIT_AVAILABLE = original
