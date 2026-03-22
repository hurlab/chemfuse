"""Tests for analyze/chemspace.py (SPEC-CF-005)."""

from __future__ import annotations

import warnings
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
    return patch("chemfuse.analyze.chemspace._RDKIT_AVAILABLE", available)


def _patch_sklearn(available=True):
    return patch("chemfuse.analyze.chemspace._SKLEARN_AVAILABLE", available)


def _patch_umap(available=True):
    return patch("chemfuse.analyze.chemspace._UMAP_AVAILABLE", available)


def _make_coords(n=5, dims=2):
    return np.random.rand(n, dims)


def _mock_fp_matrix(smiles_list, fp_type="morgan", n_bits=2048):
    """Mock fp_matrix that returns a random bit matrix."""
    return np.random.randint(0, 2, size=(len(smiles_list), n_bits), dtype=np.uint8)


class TestReduceDimensions:
    def test_requires_rdkit(self):
        """Raises OptionalDependencyError when RDKit is absent."""
        with _patch_rdkit(False):
            from chemfuse.analyze import chemspace
            with pytest.raises(OptionalDependencyError):
                chemspace.reduce_dimensions(SMILES_LIST, method="pca")

    def test_requires_sklearn(self):
        """Raises OptionalDependencyError when scikit-learn is absent."""
        with _patch_rdkit(True), _patch_sklearn(False):
            from chemfuse.analyze import chemspace
            with pytest.raises(OptionalDependencyError):
                chemspace.reduce_dimensions(SMILES_LIST, method="pca")

    def test_pca_returns_nx2_array(self):
        """PCA returns an array of shape (n, 2)."""
        mock_pca = MagicMock()
        mock_pca.fit_transform.return_value = _make_coords(5, 2)

        with _patch_rdkit(True), _patch_sklearn(True), \
             patch("chemfuse.analyze.chemspace.PCA", return_value=mock_pca), \
             patch("chemfuse.analyze.chemspace._shared_fp_matrix", side_effect=_mock_fp_matrix):

            from chemfuse.analyze import chemspace
            result = chemspace.reduce_dimensions(SMILES_LIST, method="pca")

        assert isinstance(result, np.ndarray)
        assert result.shape == (5, 2)

    def test_tsne_returns_nx2_array(self):
        """t-SNE returns an array of shape (n, 2)."""
        mock_tsne = MagicMock()
        mock_tsne.fit_transform.return_value = _make_coords(5, 2)

        with _patch_rdkit(True), _patch_sklearn(True), \
             patch("chemfuse.analyze.chemspace.TSNE", return_value=mock_tsne), \
             patch("chemfuse.analyze.chemspace._shared_fp_matrix", side_effect=_mock_fp_matrix):

            from chemfuse.analyze import chemspace
            result = chemspace.reduce_dimensions(SMILES_LIST, method="tsne")

        assert isinstance(result, np.ndarray)
        assert result.shape == (5, 2)

    def test_umap_returns_nx2_array(self):
        """UMAP returns an array of shape (n, 2) when umap-learn is available."""
        mock_umap = MagicMock()
        mock_umap.fit_transform.return_value = _make_coords(5, 2)

        with _patch_rdkit(True), _patch_sklearn(True), _patch_umap(True), \
             patch("chemfuse.analyze.chemspace.UMAP", return_value=mock_umap), \
             patch("chemfuse.analyze.chemspace._shared_fp_matrix", side_effect=_mock_fp_matrix):

            from chemfuse.analyze import chemspace
            result = chemspace.reduce_dimensions(SMILES_LIST, method="umap")

        assert isinstance(result, np.ndarray)
        assert result.shape == (5, 2)

    def test_umap_falls_back_to_pca_when_unavailable(self):
        """When umap-learn is not installed, UMAP falls back to PCA with a warning."""
        mock_pca = MagicMock()
        mock_pca.fit_transform.return_value = _make_coords(5, 2)

        with _patch_rdkit(True), _patch_sklearn(True), _patch_umap(False), \
             patch("chemfuse.analyze.chemspace.PCA", return_value=mock_pca), \
             patch("chemfuse.analyze.chemspace._shared_fp_matrix", side_effect=_mock_fp_matrix):

            from chemfuse.analyze import chemspace
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = chemspace.reduce_dimensions(SMILES_LIST, method="umap")

        assert isinstance(result, np.ndarray)
        assert result.shape[1] == 2
        assert any("umap" in str(warning.message).lower() or "PCA" in str(warning.message) for warning in w)

    def test_single_compound_returns_zeros(self):
        """Single compound returns zero array of correct shape."""
        with _patch_rdkit(True), _patch_sklearn(True), \
             patch("chemfuse.analyze.chemspace._shared_fp_matrix", side_effect=_mock_fp_matrix):

            from chemfuse.analyze import chemspace
            result = chemspace.reduce_dimensions(["CCO"], method="pca")

        assert result.shape == (1, 2)
        assert np.all(result == 0)


class TestPlotChemicalSpace:
    def test_returns_plotly_figure(self):
        """plot_chemical_space returns a plotly Figure when plotly is mocked."""
        coords = _make_coords(5, 2)

        mock_fig = MagicMock()
        mock_go = MagicMock()
        mock_go.Figure.return_value = mock_fig
        mock_go.Scatter = MagicMock()

        # Mock the plotly import inside plot_chemical_space
        with patch.dict("sys.modules", {
            "plotly": MagicMock(),
            "plotly.graph_objects": mock_go,
        }):
            from chemfuse.analyze.chemspace import plot_chemical_space
            result = plot_chemical_space(coords)

        assert result is not None

    def test_raises_when_plotly_missing(self):
        """plot_chemical_space raises OptionalDependencyError when plotly is absent."""
        import sys
        coords = _make_coords(3, 2)

        # Remove plotly from sys.modules to simulate absence
        saved_plotly = sys.modules.pop("plotly", None)
        saved_go = sys.modules.pop("plotly.graph_objects", None)
        try:
            with patch.dict("sys.modules", {"plotly": None, "plotly.graph_objects": None}):
                from chemfuse.analyze.chemspace import plot_chemical_space
                with pytest.raises((OptionalDependencyError, ImportError)):
                    plot_chemical_space(coords)
        finally:
            if saved_plotly is not None:
                sys.modules["plotly"] = saved_plotly
            if saved_go is not None:
                sys.modules["plotly.graph_objects"] = saved_go
