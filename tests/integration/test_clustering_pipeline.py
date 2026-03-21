"""Integration tests for clustering and analysis pipeline (SPEC-CF-005)."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from chemfuse.models.collection import CompoundCollection
from chemfuse.models.compound import Compound

ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
CAFFEINE_SMILES = "Cn1cnc2c1c(=O)n(c(=O)n2C)C"
IBUPROFEN_SMILES = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"

THREE_SMILES = [ASPIRIN_SMILES, CAFFEINE_SMILES, IBUPROFEN_SMILES]


def _make_compound(smiles: str, name: str, activity: float | None = None) -> Compound:
    c = Compound(smiles=smiles, name=name, sources=["test"])
    if activity is not None:
        object.__setattr__(c, "_activity", activity)
    return c


def _make_collection(*specs: tuple[str, str]) -> CompoundCollection:
    return CompoundCollection(compounds=[_make_compound(s, n) for s, n in specs])


# ---------------------------------------------------------------------------
# CompoundCollection.cluster() tests
# ---------------------------------------------------------------------------

class TestCollectionCluster:
    def test_butina_returns_labels_for_each_compound(self) -> None:
        """cluster(method='butina') returns one label per compound."""
        collection = _make_collection(
            (ASPIRIN_SMILES, "aspirin"),
            (CAFFEINE_SMILES, "caffeine"),
            (IBUPROFEN_SMILES, "ibuprofen"),
        )

        with patch("chemfuse.analyze.clustering.butina_clustering", return_value=[0, 1, 0]):
            labels = collection.cluster(method="butina", cutoff=0.4)

        assert len(labels) == 3
        assert labels == [0, 1, 0]

    def test_kmeans_returns_labels_for_each_compound(self) -> None:
        """cluster(method='kmeans') returns one label per compound."""
        collection = _make_collection(
            (ASPIRIN_SMILES, "aspirin"),
            (CAFFEINE_SMILES, "caffeine"),
            (IBUPROFEN_SMILES, "ibuprofen"),
        )

        with patch("chemfuse.analyze.clustering.kmeans_clustering", return_value=[0, 1, 2]):
            labels = collection.cluster(method="kmeans", n_clusters=3)

        assert len(labels) == 3

    def test_cluster_stores_labels_on_compounds(self) -> None:
        """cluster() attaches _cluster_label to each compound."""
        collection = _make_collection(
            (ASPIRIN_SMILES, "aspirin"),
            (CAFFEINE_SMILES, "caffeine"),
        )

        with patch("chemfuse.analyze.clustering.butina_clustering", return_value=[0, 1]):
            collection.cluster(method="butina")

        for compound, expected_label in zip(collection.compounds, [0, 1], strict=False):
            assert compound._cluster_label == expected_label

    def test_cluster_default_method_is_butina(self) -> None:
        """cluster() default method is butina."""
        collection = _make_collection(
            (ASPIRIN_SMILES, "aspirin"),
            (CAFFEINE_SMILES, "caffeine"),
        )

        with patch("chemfuse.analyze.clustering.butina_clustering", return_value=[0, 0]) as mock_bc:
            collection.cluster()

        mock_bc.assert_called_once()

    def test_cluster_passes_smiles_list(self) -> None:
        """cluster() extracts smiles from compounds and passes to algorithm."""
        collection = _make_collection(
            (ASPIRIN_SMILES, "aspirin"),
            (CAFFEINE_SMILES, "caffeine"),
        )

        with patch("chemfuse.analyze.clustering.butina_clustering", return_value=[0, 1]) as mock_bc:
            collection.cluster(method="butina", cutoff=0.5)

        call_args = mock_bc.call_args
        passed_smiles = call_args.args[0]
        assert passed_smiles == [ASPIRIN_SMILES, CAFFEINE_SMILES]

    def test_cluster_cutoff_passed_correctly(self) -> None:
        """cluster() passes cutoff parameter to butina_clustering."""
        collection = _make_collection(
            (ASPIRIN_SMILES, "aspirin"),
            (CAFFEINE_SMILES, "caffeine"),
        )

        with patch("chemfuse.analyze.clustering.butina_clustering", return_value=[0, 1]) as mock_bc:
            collection.cluster(method="butina", cutoff=0.3)

        call_kwargs = mock_bc.call_args.kwargs
        assert call_kwargs.get("cutoff") == 0.3 or (
            len(mock_bc.call_args.args) > 1 and mock_bc.call_args.args[1] == 0.3
        )

    def test_cluster_n_clusters_passed_to_kmeans(self) -> None:
        """cluster() passes n_clusters to kmeans_clustering."""
        collection = _make_collection(
            (ASPIRIN_SMILES, "aspirin"),
            (CAFFEINE_SMILES, "caffeine"),
            (IBUPROFEN_SMILES, "ibuprofen"),
        )

        with patch("chemfuse.analyze.clustering.kmeans_clustering", return_value=[0, 1, 2]) as mock_kc:
            collection.cluster(method="kmeans", n_clusters=3)

        call_kwargs = mock_kc.call_args.kwargs
        assert call_kwargs.get("n_clusters") == 3 or (
            len(mock_kc.call_args.args) > 1 and mock_kc.call_args.args[1] == 3
        )


# ---------------------------------------------------------------------------
# CompoundCollection.reduce_dimensions() tests
# ---------------------------------------------------------------------------

class TestCollectionReduceDimensions:
    def test_returns_numpy_array(self) -> None:
        """reduce_dimensions returns a numpy array."""
        collection = _make_collection(
            (ASPIRIN_SMILES, "aspirin"),
            (CAFFEINE_SMILES, "caffeine"),
            (IBUPROFEN_SMILES, "ibuprofen"),
        )
        mock_coords = np.random.rand(3, 2)

        with patch("chemfuse.analyze.chemspace.reduce_dimensions", return_value=mock_coords):
            result = collection.reduce_dimensions(method="pca")

        assert isinstance(result, np.ndarray)

    def test_returns_n_x_2_array(self) -> None:
        """reduce_dimensions returns (n, 2) array."""
        collection = _make_collection(
            (ASPIRIN_SMILES, "aspirin"),
            (CAFFEINE_SMILES, "caffeine"),
            (IBUPROFEN_SMILES, "ibuprofen"),
        )
        mock_coords = np.random.rand(3, 2)

        with patch("chemfuse.analyze.chemspace.reduce_dimensions", return_value=mock_coords):
            result = collection.reduce_dimensions(method="pca")

        assert result.shape == (3, 2)

    def test_passes_smiles_to_reduce(self) -> None:
        """reduce_dimensions passes compound SMILES list to the function."""
        collection = _make_collection(
            (ASPIRIN_SMILES, "aspirin"),
            (CAFFEINE_SMILES, "caffeine"),
        )
        mock_coords = np.random.rand(2, 2)

        with patch("chemfuse.analyze.chemspace.reduce_dimensions", return_value=mock_coords) as mock_rd:
            collection.reduce_dimensions(method="umap")

        passed_smiles = mock_rd.call_args.args[0]
        assert set(passed_smiles) == {ASPIRIN_SMILES, CAFFEINE_SMILES}


# ---------------------------------------------------------------------------
# CompoundCollection.detect_activity_cliffs() tests
# ---------------------------------------------------------------------------

class TestCollectionDetectActivityCliffs:
    def test_returns_list_of_cliffs(self) -> None:
        """detect_activity_cliffs returns a list of cliff dicts."""
        collection = _make_collection(
            (ASPIRIN_SMILES, "aspirin"),
            (CAFFEINE_SMILES, "caffeine"),
        )

        mock_cliff = {
            "smiles_1": ASPIRIN_SMILES,
            "smiles_2": CAFFEINE_SMILES,
            "similarity": 0.9,
            "activity_1": 10.0,
            "activity_2": 100.0,
            "activity_diff": 90.0,
            "cliff_score": 81.0,
        }

        with patch("chemfuse.analyze.sar.detect_activity_cliffs", return_value=[mock_cliff]):
            cliffs = collection.detect_activity_cliffs(activity_col="ic50")

        assert isinstance(cliffs, list)
        assert len(cliffs) == 1
        assert cliffs[0]["cliff_score"] == pytest.approx(81.0)

    def test_empty_collection_returns_empty(self) -> None:
        """detect_activity_cliffs on empty collection returns empty list."""
        collection = CompoundCollection(compounds=[])

        with patch("chemfuse.analyze.sar.detect_activity_cliffs", return_value=[]):
            cliffs = collection.detect_activity_cliffs()

        assert cliffs == []

    def test_passes_smiles_and_activity_col(self) -> None:
        """detect_activity_cliffs passes correct args to sar module."""
        collection = _make_collection(
            (ASPIRIN_SMILES, "aspirin"),
            (CAFFEINE_SMILES, "caffeine"),
        )

        with patch("chemfuse.analyze.sar.detect_activity_cliffs", return_value=[]) as mock_detect:
            collection.detect_activity_cliffs(activity_col="ki", sim_threshold=0.7)

        kwargs = mock_detect.call_args.kwargs
        assert kwargs.get("activity_col") == "ki"
        assert kwargs.get("sim_threshold") == pytest.approx(0.7)


# ---------------------------------------------------------------------------
# End-to-end: cluster then visualize
# ---------------------------------------------------------------------------

class TestClusterThenVisualize:
    def test_cluster_then_visualize_uses_stored_labels(self) -> None:
        """visualize_chemical_space uses _cluster_label if already clustered."""
        collection = _make_collection(
            (ASPIRIN_SMILES, "aspirin"),
            (CAFFEINE_SMILES, "caffeine"),
            (IBUPROFEN_SMILES, "ibuprofen"),
        )

        mock_coords = np.random.rand(3, 2)
        mock_fig = MagicMock()

        with patch("chemfuse.analyze.clustering.butina_clustering", return_value=[0, 1, 0]):
            collection.cluster(method="butina")

        with patch("chemfuse.analyze.chemspace.reduce_dimensions", return_value=mock_coords), \
             patch("chemfuse.analyze.chemspace.plot_chemical_space", return_value=mock_fig) as mock_plot:
            fig = collection.visualize_chemical_space(method="pca", color_by="cluster")

        assert fig is mock_fig
        # plot_chemical_space should receive the cluster labels as colors
        plot_kwargs = mock_plot.call_args.kwargs
        colors = plot_kwargs.get("colors")
        assert colors is not None
        # Colors should reflect stored cluster labels [0, 1, 0]
        assert colors == ["0", "1", "0"]
