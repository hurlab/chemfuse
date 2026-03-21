"""Tests for CompoundCollection analysis methods (SPEC-CF-005)."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from chemfuse.models.collection import CompoundCollection
from chemfuse.models.compound import Compound
from chemfuse.models.prediction import ADMETProfile

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
CAFFEINE = "Cn1cnc2c1c(=O)n(c(=O)n2C)C"
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"


def _make_compound(smiles: str, name: str) -> Compound:
    return Compound(smiles=smiles, name=name, sources=["test"])


def _make_collection(*pairs: tuple[str, str]) -> CompoundCollection:
    return CompoundCollection(compounds=[_make_compound(s, n) for s, n in pairs])


# ---------------------------------------------------------------------------
# predict_admet method
# ---------------------------------------------------------------------------

class TestPredictAdmetMethod:
    def test_method_exists(self) -> None:
        """CompoundCollection has predict_admet method."""
        collection = _make_collection((ASPIRIN, "aspirin"))
        assert hasattr(collection, "predict_admet")
        assert callable(collection.predict_admet)

    def test_returns_none(self) -> None:
        """predict_admet returns None (in-place operation)."""
        collection = _make_collection((ASPIRIN, "aspirin"))
        with patch("chemfuse.compute.admet.predict_admet", return_value=ADMETProfile(smiles=ASPIRIN)):
            result = collection.predict_admet()
        assert result is None

    def test_stores_profile_attribute(self) -> None:
        """predict_admet stores _admet_profile on each compound."""
        collection = _make_collection((ASPIRIN, "aspirin"))
        profile = ADMETProfile(smiles=ASPIRIN, overall_score=0.8)

        with patch("chemfuse.compute.admet.predict_admet", return_value=profile):
            collection.predict_admet()

        compound = collection.compounds[0]
        assert hasattr(compound, "_admet_profile")
        stored = compound._admet_profile
        assert stored is profile

    def test_called_once_per_compound(self) -> None:
        """predict_admet calls the prediction function once per compound."""
        collection = _make_collection(
            (ASPIRIN, "aspirin"),
            (CAFFEINE, "caffeine"),
            (IBUPROFEN, "ibuprofen"),
        )
        profile = ADMETProfile(smiles=ASPIRIN)

        with patch("chemfuse.compute.admet.predict_admet", return_value=profile) as mock_pred:
            collection.predict_admet()

        assert mock_pred.call_count == 3

    def test_empty_collection_does_nothing(self) -> None:
        """predict_admet on empty collection makes no prediction calls."""
        collection = CompoundCollection(compounds=[])
        with patch("chemfuse.compute.admet.predict_admet") as mock_pred:
            collection.predict_admet()
        mock_pred.assert_not_called()


# ---------------------------------------------------------------------------
# cluster method
# ---------------------------------------------------------------------------

class TestClusterMethod:
    def test_method_exists(self) -> None:
        """CompoundCollection has cluster method."""
        collection = _make_collection((ASPIRIN, "aspirin"))
        assert hasattr(collection, "cluster")
        assert callable(collection.cluster)

    def test_returns_list_of_ints(self) -> None:
        """cluster() returns list of integer labels."""
        collection = _make_collection((ASPIRIN, "aspirin"), (CAFFEINE, "caffeine"))
        with patch("chemfuse.analyze.clustering.butina_clustering", return_value=[0, 1]):
            labels = collection.cluster()
        assert isinstance(labels, list)
        assert all(isinstance(lbl, int) for lbl in labels)

    def test_length_matches_collection(self) -> None:
        """cluster() returns same number of labels as compounds."""
        collection = _make_collection(
            (ASPIRIN, "aspirin"),
            (CAFFEINE, "caffeine"),
            (IBUPROFEN, "ibuprofen"),
        )
        with patch("chemfuse.analyze.clustering.butina_clustering", return_value=[0, 1, 0]):
            labels = collection.cluster()
        assert len(labels) == 3

    def test_butina_method(self) -> None:
        """cluster(method='butina') delegates to butina_clustering."""
        collection = _make_collection((ASPIRIN, "aspirin"), (CAFFEINE, "caffeine"))
        with patch("chemfuse.analyze.clustering.butina_clustering", return_value=[0, 0]) as mock_bc:
            collection.cluster(method="butina")
        mock_bc.assert_called_once()

    def test_kmeans_method(self) -> None:
        """cluster(method='kmeans') delegates to kmeans_clustering."""
        collection = _make_collection((ASPIRIN, "aspirin"), (CAFFEINE, "caffeine"))
        with patch("chemfuse.analyze.clustering.kmeans_clustering", return_value=[0, 1]) as mock_kc:
            collection.cluster(method="kmeans")
        mock_kc.assert_called_once()

    def test_labels_attached_to_compounds(self) -> None:
        """cluster() stores _cluster_label on each compound."""
        collection = _make_collection(
            (ASPIRIN, "aspirin"),
            (CAFFEINE, "caffeine"),
        )
        with patch("chemfuse.analyze.clustering.butina_clustering", return_value=[0, 1]):
            collection.cluster(method="butina")
        assert collection.compounds[0]._cluster_label == 0
        assert collection.compounds[1]._cluster_label == 1


# ---------------------------------------------------------------------------
# reduce_dimensions method
# ---------------------------------------------------------------------------

class TestReduceDimensionsMethod:
    def test_method_exists(self) -> None:
        """CompoundCollection has reduce_dimensions method."""
        collection = _make_collection((ASPIRIN, "aspirin"))
        assert hasattr(collection, "reduce_dimensions")

    def test_returns_ndarray(self) -> None:
        """reduce_dimensions returns a numpy array."""
        collection = _make_collection((ASPIRIN, "aspirin"), (CAFFEINE, "caffeine"))
        mock_coords = np.zeros((2, 2))
        with patch("chemfuse.analyze.chemspace.reduce_dimensions", return_value=mock_coords):
            result = collection.reduce_dimensions(method="pca")
        assert isinstance(result, np.ndarray)

    def test_array_shape_matches_collection(self) -> None:
        """reduce_dimensions array has same row count as compounds."""
        collection = _make_collection(
            (ASPIRIN, "aspirin"),
            (CAFFEINE, "caffeine"),
            (IBUPROFEN, "ibuprofen"),
        )
        mock_coords = np.zeros((3, 2))
        with patch("chemfuse.analyze.chemspace.reduce_dimensions", return_value=mock_coords):
            result = collection.reduce_dimensions()
        assert result.shape[0] == 3


# ---------------------------------------------------------------------------
# detect_activity_cliffs method
# ---------------------------------------------------------------------------

class TestDetectActivityCliffsMethod:
    def test_method_exists(self) -> None:
        """CompoundCollection has detect_activity_cliffs method."""
        collection = _make_collection((ASPIRIN, "aspirin"))
        assert hasattr(collection, "detect_activity_cliffs")

    def test_returns_list(self) -> None:
        """detect_activity_cliffs returns a list."""
        collection = _make_collection(
            (ASPIRIN, "aspirin"),
            (CAFFEINE, "caffeine"),
        )
        with patch("chemfuse.analyze.sar.detect_activity_cliffs", return_value=[]):
            result = collection.detect_activity_cliffs()
        assert isinstance(result, list)

    def test_passes_activity_col(self) -> None:
        """detect_activity_cliffs passes activity_col argument."""
        collection = _make_collection(
            (ASPIRIN, "aspirin"),
            (CAFFEINE, "caffeine"),
        )
        with patch("chemfuse.analyze.sar.detect_activity_cliffs", return_value=[]) as mock_detect:
            collection.detect_activity_cliffs(activity_col="ki")
        kwargs = mock_detect.call_args.kwargs
        assert kwargs.get("activity_col") == "ki"

    def test_passes_sim_threshold(self) -> None:
        """detect_activity_cliffs passes sim_threshold argument."""
        collection = _make_collection(
            (ASPIRIN, "aspirin"),
            (CAFFEINE, "caffeine"),
        )
        with patch("chemfuse.analyze.sar.detect_activity_cliffs", return_value=[]) as mock_detect:
            collection.detect_activity_cliffs(sim_threshold=0.6)
        kwargs = mock_detect.call_args.kwargs
        assert kwargs.get("sim_threshold") == pytest.approx(0.6)

    def test_compound_smiles_passed_in_dicts(self) -> None:
        """detect_activity_cliffs passes each compound's smiles in the dict list."""
        collection = _make_collection(
            (ASPIRIN, "aspirin"),
            (CAFFEINE, "caffeine"),
        )
        with patch("chemfuse.analyze.sar.detect_activity_cliffs", return_value=[]) as mock_detect:
            collection.detect_activity_cliffs()

        compound_dicts = mock_detect.call_args.args[0]
        smiles_in_dicts = [d["smiles"] for d in compound_dicts]
        assert ASPIRIN in smiles_in_dicts
        assert CAFFEINE in smiles_in_dicts


# ---------------------------------------------------------------------------
# visualize_chemical_space method
# ---------------------------------------------------------------------------

class TestVisualizeChemicalSpace:
    def test_method_exists(self) -> None:
        """CompoundCollection has visualize_chemical_space method."""
        collection = _make_collection((ASPIRIN, "aspirin"))
        assert hasattr(collection, "visualize_chemical_space")

    def test_returns_figure(self) -> None:
        """visualize_chemical_space returns a plotly-like figure object."""
        collection = _make_collection(
            (ASPIRIN, "aspirin"),
            (CAFFEINE, "caffeine"),
        )
        mock_coords = np.zeros((2, 2))
        mock_fig = MagicMock()

        with patch("chemfuse.analyze.chemspace.reduce_dimensions", return_value=mock_coords), \
             patch("chemfuse.analyze.chemspace.plot_chemical_space", return_value=mock_fig):
            fig = collection.visualize_chemical_space(method="pca")

        assert fig is mock_fig

    def test_color_by_cluster_uses_stored_labels(self) -> None:
        """visualize_chemical_space uses _cluster_label when color_by='cluster'."""
        collection = _make_collection(
            (ASPIRIN, "aspirin"),
            (CAFFEINE, "caffeine"),
        )
        # Pre-attach cluster labels
        object.__setattr__(collection.compounds[0], "_cluster_label", 0)
        object.__setattr__(collection.compounds[1], "_cluster_label", 1)

        mock_coords = np.zeros((2, 2))
        mock_fig = MagicMock()

        with patch("chemfuse.analyze.chemspace.reduce_dimensions", return_value=mock_coords), \
             patch("chemfuse.analyze.chemspace.plot_chemical_space", return_value=mock_fig) as mock_plot:
            collection.visualize_chemical_space(color_by="cluster")

        kwargs = mock_plot.call_args.kwargs
        colors = kwargs.get("colors")
        assert colors == ["0", "1"]

    def test_color_by_cluster_no_labels_passes_none_colors(self) -> None:
        """visualize_chemical_space passes None colors when no cluster labels set."""
        collection = _make_collection(
            (ASPIRIN, "aspirin"),
            (CAFFEINE, "caffeine"),
        )

        mock_coords = np.zeros((2, 2))
        mock_fig = MagicMock()

        with patch("chemfuse.analyze.chemspace.reduce_dimensions", return_value=mock_coords), \
             patch("chemfuse.analyze.chemspace.plot_chemical_space", return_value=mock_fig) as mock_plot:
            collection.visualize_chemical_space(color_by="cluster")

        kwargs = mock_plot.call_args.kwargs
        colors = kwargs.get("colors")
        assert colors is None
