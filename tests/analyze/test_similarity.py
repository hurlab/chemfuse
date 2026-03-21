"""Tests for analyze/similarity.py (SPEC-CF-005)."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from chemfuse.core.exceptions import OptionalDependencyError

SMILES_A = "CC(=O)Oc1ccccc1C(=O)O"  # aspirin
SMILES_B = "Cn1cnc2c1c(=O)n(c(=O)n2C)C"  # caffeine


def _make_fp_mock(bits=2048):
    fp = MagicMock()
    return fp


def _patch_rdkit(available=True):
    """Return a patch context for rdkit availability in the similarity module."""
    return patch("chemfuse.analyze.similarity._RDKIT_AVAILABLE", available)


class TestTanimotoSimilarity:
    def test_requires_rdkit(self):
        """Raises OptionalDependencyError when RDKit is absent."""
        with _patch_rdkit(False):
            from chemfuse.analyze import similarity
            with pytest.raises(OptionalDependencyError):
                similarity.tanimoto_similarity(SMILES_A, SMILES_B)

    def test_identical_smiles_returns_one(self):
        """tanimoto_similarity for identical SMILES returns 1.0."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem, \
             patch("chemfuse.analyze.similarity.DataStructs") as mock_ds, \
             patch("chemfuse.analyze.similarity.AllChem") as mock_allchem:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            fp = _make_fp_mock()
            mock_allchem.GetMorganFingerprintAsBitVect.return_value = fp
            mock_ds.TanimotoSimilarity.return_value = 1.0

            from chemfuse.analyze import similarity
            result = similarity.tanimoto_similarity(SMILES_A, SMILES_A)
            assert result == 1.0

    def test_different_smiles_returns_float(self):
        """tanimoto_similarity returns a float in [0, 1]."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem, \
             patch("chemfuse.analyze.similarity.DataStructs") as mock_ds, \
             patch("chemfuse.analyze.similarity.AllChem") as mock_allchem:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            fp = _make_fp_mock()
            mock_allchem.GetMorganFingerprintAsBitVect.return_value = fp
            mock_ds.TanimotoSimilarity.return_value = 0.35

            from chemfuse.analyze import similarity
            result = similarity.tanimoto_similarity(SMILES_A, SMILES_B)
            assert isinstance(result, float)
            assert 0.0 <= result <= 1.0

    def test_invalid_smiles_raises_value_error(self):
        """tanimoto_similarity raises ValueError for invalid SMILES."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem:

            mock_chem.MolFromSmiles.return_value = None

            from chemfuse.analyze import similarity
            with pytest.raises(ValueError, match="Invalid SMILES"):
                similarity.tanimoto_similarity("BADSMILES", SMILES_B)


class TestBulkTanimoto:
    def test_requires_rdkit(self):
        """Raises OptionalDependencyError when RDKit is absent."""
        with _patch_rdkit(False):
            from chemfuse.analyze import similarity
            with pytest.raises(OptionalDependencyError):
                similarity.bulk_tanimoto(SMILES_A, [SMILES_B])

    def test_returns_sorted_descending(self):
        """bulk_tanimoto returns results sorted by descending similarity."""
        sims = [0.5, 0.9, 0.1]

        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem, \
             patch("chemfuse.analyze.similarity.DataStructs") as mock_ds, \
             patch("chemfuse.analyze.similarity.AllChem") as mock_allchem:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            fp = _make_fp_mock()
            mock_allchem.GetMorganFingerprintAsBitVect.return_value = fp
            mock_ds.TanimotoSimilarity.side_effect = sims

            from chemfuse.analyze import similarity
            results = similarity.bulk_tanimoto(SMILES_A, [SMILES_A, SMILES_B, "CCO"])

        assert len(results) == 3
        scores = [r[1] for r in results]
        assert scores == sorted(scores, reverse=True)

    def test_empty_targets_returns_empty(self):
        """bulk_tanimoto with empty targets returns empty list."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem, \
             patch("chemfuse.analyze.similarity.AllChem") as mock_allchem:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_allchem.GetMorganFingerprintAsBitVect.return_value = _make_fp_mock()

            from chemfuse.analyze import similarity
            results = similarity.bulk_tanimoto(SMILES_A, [])

        assert results == []


class TestTanimotoMatrix:
    def test_requires_rdkit(self):
        """Raises OptionalDependencyError when RDKit is absent."""
        with _patch_rdkit(False):
            from chemfuse.analyze import similarity
            with pytest.raises(OptionalDependencyError):
                similarity.tanimoto_matrix([SMILES_A])

    def test_matrix_shape(self):
        """tanimoto_matrix returns n x n array."""
        smiles_list = [SMILES_A, SMILES_B, "CCO"]
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem, \
             patch("chemfuse.analyze.similarity.DataStructs") as mock_ds, \
             patch("chemfuse.analyze.similarity.AllChem") as mock_allchem:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_allchem.GetMorganFingerprintAsBitVect.return_value = _make_fp_mock()
            mock_ds.TanimotoSimilarity.return_value = 0.4

            from chemfuse.analyze import similarity
            mat = similarity.tanimoto_matrix(smiles_list)

        assert isinstance(mat, np.ndarray)
        assert mat.shape == (3, 3)

    def test_matrix_diagonal_is_one(self):
        """Diagonal of tanimoto_matrix should be 1.0 (self-similarity)."""
        smiles_list = [SMILES_A, SMILES_B]

        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem, \
             patch("chemfuse.analyze.similarity.DataStructs") as mock_ds, \
             patch("chemfuse.analyze.similarity.AllChem") as mock_allchem:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_allchem.GetMorganFingerprintAsBitVect.return_value = _make_fp_mock()
            mock_ds.TanimotoSimilarity.return_value = 0.35

            from chemfuse.analyze import similarity
            mat = similarity.tanimoto_matrix(smiles_list)

        assert mat[0, 0] == pytest.approx(1.0)
        assert mat[1, 1] == pytest.approx(1.0)

    def test_matrix_is_symmetric(self):
        """tanimoto_matrix must be symmetric."""
        smiles_list = [SMILES_A, SMILES_B, "CCO"]

        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem, \
             patch("chemfuse.analyze.similarity.DataStructs") as mock_ds, \
             patch("chemfuse.analyze.similarity.AllChem") as mock_allchem:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_allchem.GetMorganFingerprintAsBitVect.return_value = _make_fp_mock()
            mock_ds.TanimotoSimilarity.return_value = 0.4

            from chemfuse.analyze import similarity
            mat = similarity.tanimoto_matrix(smiles_list)

        np.testing.assert_array_equal(mat, mat.T)


class TestFindNearestNeighbors:
    def test_returns_top_n(self):
        """find_nearest_neighbors returns at most n results."""
        targets = [SMILES_A, SMILES_B, "CCO", "c1ccccc1", "CCC"]

        with _patch_rdkit(True), \
             patch("chemfuse.analyze.similarity.Chem") as mock_chem, \
             patch("chemfuse.analyze.similarity.DataStructs") as mock_ds, \
             patch("chemfuse.analyze.similarity.AllChem") as mock_allchem:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_allchem.GetMorganFingerprintAsBitVect.return_value = _make_fp_mock()
            mock_ds.TanimotoSimilarity.return_value = 0.5

            from chemfuse.analyze import similarity
            results = similarity.find_nearest_neighbors(SMILES_A, targets, n=3)

        assert len(results) == 3
