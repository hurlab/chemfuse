"""Tests for analyze/clustering.py (SPEC-CF-005)."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from chemfuse.core.exceptions import OptionalDependencyError

SMILES_LIST = [
    "CC(=O)Oc1ccccc1C(=O)O",
    "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "CCO",
    "c1ccccc1",
    "CCC",
]


def _patch_rdkit(available=True):
    return patch("chemfuse.analyze.clustering._RDKIT_AVAILABLE", available)


def _patch_sklearn(available=True):
    return patch("chemfuse.analyze.clustering._SKLEARN_AVAILABLE", available)


def _make_fp_mock():
    return MagicMock()


class TestButinaClustering:
    def test_requires_rdkit(self):
        """Raises OptionalDependencyError when RDKit is absent."""
        with _patch_rdkit(False):
            from chemfuse.analyze import clustering
            with pytest.raises(OptionalDependencyError):
                clustering.butina_clustering(SMILES_LIST)

    def test_returns_label_per_smiles(self):
        """butina_clustering returns one label per input SMILES."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.clustering.Chem") as mock_chem, \
             patch("chemfuse.analyze.clustering.DataStructs") as mock_ds, \
             patch("chemfuse.analyze.clustering.AllChem") as mock_allchem, \
             patch("chemfuse.analyze.clustering.Butina") as mock_butina:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_allchem.GetMorganFingerprintAsBitVect.return_value = _make_fp_mock()
            mock_ds.TanimotoSimilarity.return_value = 0.3

            # Return clusters: cluster 0 has members [0,1,2,3,4]
            mock_butina.ClusterData.return_value = [(0, 1, 2), (3,), (4,)]

            from chemfuse.analyze import clustering
            labels = clustering.butina_clustering(SMILES_LIST)

        assert len(labels) == len(SMILES_LIST)

    def test_returns_integer_labels(self):
        """butina_clustering returns integer labels."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.clustering.Chem") as mock_chem, \
             patch("chemfuse.analyze.clustering.DataStructs") as mock_ds, \
             patch("chemfuse.analyze.clustering.AllChem") as mock_allchem, \
             patch("chemfuse.analyze.clustering.Butina") as mock_butina:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_allchem.GetMorganFingerprintAsBitVect.return_value = _make_fp_mock()
            mock_ds.TanimotoSimilarity.return_value = 0.2
            mock_butina.ClusterData.return_value = [(0, 1, 2, 3, 4),]

            from chemfuse.analyze import clustering
            labels = clustering.butina_clustering(SMILES_LIST)

        assert all(isinstance(lbl, int) for lbl in labels)

    def test_empty_input_returns_empty(self):
        """butina_clustering with empty input returns empty list."""
        with _patch_rdkit(True):
            from chemfuse.analyze import clustering
            labels = clustering.butina_clustering([])
        assert labels == []

    def test_single_compound(self):
        """butina_clustering with single compound returns [0]."""
        with _patch_rdkit(True), \
             patch("chemfuse.analyze.clustering.Chem") as mock_chem, \
             patch("chemfuse.analyze.clustering.AllChem") as mock_allchem:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_allchem.GetMorganFingerprintAsBitVect.return_value = _make_fp_mock()

            from chemfuse.analyze import clustering
            labels = clustering.butina_clustering(["CCO"])

        assert len(labels) == 1


class TestKMeansClustering:
    def test_requires_rdkit(self):
        """Raises OptionalDependencyError when RDKit is absent."""
        with _patch_rdkit(False):
            from chemfuse.analyze import clustering
            with pytest.raises(OptionalDependencyError):
                clustering.kmeans_clustering(SMILES_LIST, n_clusters=2)

    def test_requires_sklearn(self):
        """Raises OptionalDependencyError when scikit-learn is absent."""
        with _patch_rdkit(True), _patch_sklearn(False):
            from chemfuse.analyze import clustering
            with pytest.raises(OptionalDependencyError):
                clustering.kmeans_clustering(SMILES_LIST, n_clusters=2)

    def test_returns_label_per_smiles(self):
        """kmeans_clustering returns one label per input SMILES."""
        mock_km = MagicMock()
        mock_km.fit_predict.return_value = np.array([0, 1, 2, 0, 1])

        with _patch_rdkit(True), _patch_sklearn(True), \
             patch("chemfuse.analyze.clustering.Chem") as mock_chem, \
             patch("chemfuse.analyze.clustering.DataStructs") as mock_ds, \
             patch("chemfuse.analyze.clustering.AllChem") as mock_allchem, \
             patch("chemfuse.analyze.clustering.MACCSkeys"), \
             patch("chemfuse.analyze.clustering.KMeans", return_value=mock_km):

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            fp = MagicMock()
            mock_allchem.GetMorganFingerprintAsBitVect.return_value = fp
            mock_ds.ConvertToNumpyArray = lambda fp, arr: None

            from chemfuse.analyze import clustering
            labels = clustering.kmeans_clustering(SMILES_LIST, n_clusters=3)

        assert len(labels) == len(SMILES_LIST)

    def test_returns_exactly_k_unique_labels(self):
        """kmeans_clustering with n_clusters=5 returns exactly 5 unique labels."""
        mock_km = MagicMock()
        mock_km.fit_predict.return_value = np.array([0, 1, 2, 3, 4])

        with _patch_rdkit(True), _patch_sklearn(True), \
             patch("chemfuse.analyze.clustering.Chem") as mock_chem, \
             patch("chemfuse.analyze.clustering.DataStructs") as mock_ds, \
             patch("chemfuse.analyze.clustering.AllChem") as mock_allchem, \
             patch("chemfuse.analyze.clustering.KMeans", return_value=mock_km):

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            fp = MagicMock()
            mock_allchem.GetMorganFingerprintAsBitVect.return_value = fp
            mock_ds.ConvertToNumpyArray = lambda fp, arr: None

            from chemfuse.analyze import clustering
            labels = clustering.kmeans_clustering(SMILES_LIST * 10, n_clusters=5)

        assert len(set(labels)) == 5

    def test_empty_input(self):
        """kmeans_clustering with empty input returns empty list."""
        with _patch_rdkit(True), _patch_sklearn(True):
            from chemfuse.analyze import clustering
            assert clustering.kmeans_clustering([]) == []


class TestSilhouetteScore:
    def test_returns_none_without_sklearn(self):
        """compute_silhouette returns None when sklearn is absent."""
        with _patch_sklearn(False):
            from chemfuse.analyze import clustering
            result = clustering.compute_silhouette(SMILES_LIST, [0, 0, 1, 1, 0])
        assert result is None

    def test_returns_none_with_single_cluster(self):
        """compute_silhouette returns None when all labels are the same."""
        with _patch_rdkit(True), _patch_sklearn(True):
            from chemfuse.analyze import clustering
            result = clustering.compute_silhouette(SMILES_LIST, [0, 0, 0, 0, 0])
        assert result is None

    def test_returns_float_with_valid_input(self):
        """compute_silhouette returns a float score with valid clustering."""
        from unittest.mock import patch

        with _patch_rdkit(True), _patch_sklearn(True), \
             patch("chemfuse.analyze.clustering.silhouette_score", return_value=0.42), \
             patch("chemfuse.analyze.clustering.Chem") as mock_chem, \
             patch("chemfuse.analyze.clustering.DataStructs") as mock_ds, \
             patch("chemfuse.analyze.clustering.AllChem") as mock_allchem:

            mol = MagicMock()
            mock_chem.MolFromSmiles.return_value = mol
            mock_allchem.GetMorganFingerprintAsBitVect.return_value = MagicMock()
            mock_ds.ConvertToNumpyArray = lambda fp, arr: None

            from chemfuse.analyze import clustering
            result = clustering.compute_silhouette(SMILES_LIST, [0, 0, 1, 1, 0])

        assert isinstance(result, float)
        assert result == pytest.approx(0.42)
