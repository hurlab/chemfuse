"""Tests for analyze/landscape.py (CF-E10: SAR Landscape Visualization)."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import numpy as np
import pytest

SMILES = ["CCO", "CCCO", "CCCCO", "c1ccccc1", "c1ccncc1", "CC(=O)O"]
ACTIVITIES = [5.0, 5.5, 6.0, 7.5, 8.0, 4.5]
_N = len(SMILES)


def _mock_reduce_dimensions(smiles_list, method="pca", fp_type="morgan", n_components=2, **kwargs):
    """Mock reduce_dimensions returning deterministic zero coords."""
    n = len(smiles_list)
    if n == 0:
        return np.empty((0, 2), dtype=float)
    return np.zeros((n, 2), dtype=float)


def _mock_tanimoto_matrix(smiles_list, fp_type="morgan"):
    """Mock tanimoto_matrix returning identity (high similarity within alcohols)."""
    n = len(smiles_list)
    mat = np.eye(n, dtype=float)
    # Make adjacent pairs highly similar (alcohols / chains)
    for i in range(n - 1):
        mat[i, i + 1] = 0.9
        mat[i + 1, i] = 0.9
    return mat


class TestSarLandscape:
    """Tests for the sar_landscape() function."""

    def test_returns_expected_keys(self):
        """sar_landscape returns dict with all required keys."""
        from chemfuse.analyze.landscape import sar_landscape

        with (
            patch("chemfuse.analyze.landscape.reduce_dimensions", side_effect=_mock_reduce_dimensions),
            patch("chemfuse.analyze.landscape._detect_cliff_pairs", return_value=[]),
            patch("chemfuse.analyze.landscape._build_plotly_figure", return_value=None),
        ):
            result = sar_landscape(SMILES, ACTIVITIES, method="pca")

        assert set(result.keys()) == {"coords", "activities", "smiles", "cliff_pairs", "figure"}

    def test_coords_shape(self):
        """coords is a numpy array of shape (n, 2)."""
        from chemfuse.analyze.landscape import sar_landscape

        with (
            patch("chemfuse.analyze.landscape.reduce_dimensions", side_effect=_mock_reduce_dimensions),
            patch("chemfuse.analyze.landscape._detect_cliff_pairs", return_value=[]),
            patch("chemfuse.analyze.landscape._build_plotly_figure", return_value=None),
        ):
            result = sar_landscape(SMILES, ACTIVITIES, method="pca")

        assert isinstance(result["coords"], np.ndarray)
        assert result["coords"].shape == (_N, 2)

    def test_activities_and_smiles_preserved(self):
        """activities and smiles in result match inputs."""
        from chemfuse.analyze.landscape import sar_landscape

        with (
            patch("chemfuse.analyze.landscape.reduce_dimensions", side_effect=_mock_reduce_dimensions),
            patch("chemfuse.analyze.landscape._detect_cliff_pairs", return_value=[]),
            patch("chemfuse.analyze.landscape._build_plotly_figure", return_value=None),
        ):
            result = sar_landscape(SMILES, ACTIVITIES, method="pca")

        assert result["activities"] == list(ACTIVITIES)
        assert result["smiles"] == list(SMILES)

    def test_cliff_pairs_are_list_of_dicts(self):
        """cliff_pairs is a list of dicts with expected keys."""
        from chemfuse.analyze.landscape import sar_landscape

        fake_pairs = [{"i": 0, "j": 1, "similarity": 0.9, "activity_diff": 0.5}]

        with (
            patch("chemfuse.analyze.landscape.reduce_dimensions", side_effect=_mock_reduce_dimensions),
            patch("chemfuse.analyze.landscape._detect_cliff_pairs", return_value=fake_pairs),
            patch("chemfuse.analyze.landscape._build_plotly_figure", return_value=None),
        ):
            result = sar_landscape(SMILES, ACTIVITIES, method="pca")

        assert isinstance(result["cliff_pairs"], list)
        for pair in result["cliff_pairs"]:
            assert "i" in pair
            assert "j" in pair
            assert "similarity" in pair
            assert "activity_diff" in pair

    def test_cliff_pairs_respect_threshold(self):
        """cliff_pairs only contains pairs above both thresholds (integration with _detect_cliff_pairs)."""
        from chemfuse.analyze.landscape import _detect_cliff_pairs

        with patch("chemfuse.analyze.landscape.tanimoto_matrix", side_effect=_mock_tanimoto_matrix):
            # Very high threshold: no pairs should qualify
            pairs = _detect_cliff_pairs(
                SMILES,
                ACTIVITIES,
                fp_type="morgan",
                cliff_threshold=0.999,
                activity_diff_threshold=100.0,
            )
        assert pairs == []

    def test_cliff_pairs_indices_are_valid(self):
        """All cliff pair indices are within valid range."""
        from chemfuse.analyze.landscape import _detect_cliff_pairs

        with patch("chemfuse.analyze.landscape.tanimoto_matrix", side_effect=_mock_tanimoto_matrix):
            pairs = _detect_cliff_pairs(
                SMILES,
                ACTIVITIES,
                fp_type="morgan",
                cliff_threshold=0.0,
                activity_diff_threshold=0.1,
            )

        n = len(SMILES)
        for pair in pairs:
            assert 0 <= pair["i"] < n
            assert 0 <= pair["j"] < n
            assert pair["i"] < pair["j"]

    def test_pca_method_called(self):
        """reduce_dimensions is called with method='pca'."""
        from chemfuse.analyze.landscape import sar_landscape

        mock_reduce = MagicMock(side_effect=_mock_reduce_dimensions)

        with (
            patch("chemfuse.analyze.landscape.reduce_dimensions", mock_reduce),
            patch("chemfuse.analyze.landscape._detect_cliff_pairs", return_value=[]),
            patch("chemfuse.analyze.landscape._build_plotly_figure", return_value=None),
        ):
            sar_landscape(SMILES, ACTIVITIES, method="pca")

        mock_reduce.assert_called_once()
        _, kwargs = mock_reduce.call_args
        # method may be positional or keyword
        call_args = mock_reduce.call_args
        assert "pca" in str(call_args)

    def test_empty_list_returns_empty_result(self):
        """Empty SMILES list returns a valid empty result without calling reduce_dimensions."""
        from chemfuse.analyze.landscape import sar_landscape

        result = sar_landscape([], [], method="pca")

        assert result["coords"].shape == (0, 2)
        assert result["activities"] == []
        assert result["smiles"] == []
        assert result["cliff_pairs"] == []
        assert result["figure"] is None

    def test_mismatched_lengths_raises_value_error(self):
        """Mismatched smiles_list and activities lengths raises ValueError."""
        from chemfuse.analyze.landscape import sar_landscape

        with pytest.raises(ValueError, match="same length"):
            sar_landscape(["CCO", "CCCO"], [5.0], method="pca")

    def test_figure_is_none_when_plotly_not_installed(self):
        """figure is None when plotly is not available."""
        from chemfuse.analyze.landscape import _build_plotly_figure

        coords = np.zeros((3, 2))
        with patch.dict("sys.modules", {"plotly": None, "plotly.graph_objects": None}):
            result = _build_plotly_figure(
                coords=coords,
                activities=[1.0, 2.0, 3.0],
                smiles_list=["CCO", "CCCO", "CCCCO"],
                cliff_pairs=[],
                colormap="RdYlGn_r",
                method="pca",
            )
        assert result is None

    def test_figure_in_result_when_build_returns_mock(self):
        """figure in result comes from _build_plotly_figure return value."""
        from chemfuse.analyze.landscape import sar_landscape

        fake_fig = MagicMock(name="plotly.Figure")

        with (
            patch("chemfuse.analyze.landscape.reduce_dimensions", side_effect=_mock_reduce_dimensions),
            patch("chemfuse.analyze.landscape._detect_cliff_pairs", return_value=[]),
            patch("chemfuse.analyze.landscape._build_plotly_figure", return_value=fake_fig),
        ):
            result = sar_landscape(SMILES, ACTIVITIES, method="pca")

        assert result["figure"] is fake_fig

    def test_single_compound_returns_one_row_coords(self):
        """Single compound returns shape (1, 2) coords and empty cliff_pairs."""
        from chemfuse.analyze.landscape import sar_landscape

        with (
            patch("chemfuse.analyze.landscape.reduce_dimensions", side_effect=_mock_reduce_dimensions),
            patch("chemfuse.analyze.landscape._build_plotly_figure", return_value=None),
        ):
            result = sar_landscape(["CCO"], [5.0], method="pca")

        assert result["coords"].shape == (1, 2)
        assert result["cliff_pairs"] == []

    def test_cliff_pairs_activity_diff_meets_threshold(self):
        """All returned cliff pairs meet the activity_diff_threshold."""
        from chemfuse.analyze.landscape import _detect_cliff_pairs

        activity_diff_threshold = 0.4
        with patch("chemfuse.analyze.landscape.tanimoto_matrix", side_effect=_mock_tanimoto_matrix):
            pairs = _detect_cliff_pairs(
                SMILES,
                ACTIVITIES,
                fp_type="morgan",
                cliff_threshold=0.0,
                activity_diff_threshold=activity_diff_threshold,
            )

        for pair in pairs:
            assert pair["activity_diff"] >= activity_diff_threshold


class TestDetectCliffPairs:
    """Direct unit tests for _detect_cliff_pairs helper."""

    def test_empty_returns_empty(self):
        """Empty input returns empty list."""
        from chemfuse.analyze.landscape import _detect_cliff_pairs

        with patch("chemfuse.analyze.landscape.tanimoto_matrix", side_effect=_mock_tanimoto_matrix):
            result = _detect_cliff_pairs([], [], fp_type="morgan", cliff_threshold=0.8, activity_diff_threshold=1.0)
        assert result == []

    def test_single_compound_returns_empty(self):
        """Single compound cannot form pairs."""
        from chemfuse.analyze.landscape import _detect_cliff_pairs

        with patch("chemfuse.analyze.landscape.tanimoto_matrix", side_effect=_mock_tanimoto_matrix):
            result = _detect_cliff_pairs(["CCO"], [5.0], fp_type="morgan", cliff_threshold=0.8, activity_diff_threshold=1.0)
        assert result == []

    def test_high_sim_low_diff_excluded(self):
        """Pairs with high similarity but low activity diff are excluded."""
        from chemfuse.analyze.landscape import _detect_cliff_pairs

        def mock_mat(smiles_list, fp_type="morgan"):
            return np.array([[1.0, 0.95], [0.95, 1.0]])

        # Same activity => diff = 0.0 which is < 1.0 threshold
        with patch("chemfuse.analyze.landscape.tanimoto_matrix", mock_mat):
            result = _detect_cliff_pairs(
                ["CCO", "CCCO"],
                [5.0, 5.0],
                fp_type="morgan",
                cliff_threshold=0.8,
                activity_diff_threshold=1.0,
            )
        assert result == []

    def test_high_sim_high_diff_included(self):
        """Pairs meeting both thresholds are included."""
        from chemfuse.analyze.landscape import _detect_cliff_pairs

        def mock_mat(smiles_list, fp_type="morgan"):
            return np.array([[1.0, 0.95], [0.95, 1.0]])

        with patch("chemfuse.analyze.landscape.tanimoto_matrix", mock_mat):
            result = _detect_cliff_pairs(
                ["CCO", "CCCO"],
                [5.0, 8.0],
                fp_type="morgan",
                cliff_threshold=0.8,
                activity_diff_threshold=1.0,
            )

        assert len(result) == 1
        assert result[0]["i"] == 0
        assert result[0]["j"] == 1
        assert result[0]["similarity"] == 0.95
        assert result[0]["activity_diff"] == 3.0


class TestCompoundCollectionSarLandscape:
    """Tests for CompoundCollection.plot_sar_landscape()."""

    def _make_collection_with_bioactivities(self):
        """Build a minimal CompoundCollection with bioactivity data."""
        from chemfuse.models.bioactivity import Bioactivity
        from chemfuse.models.collection import CompoundCollection
        from chemfuse.models.compound import Compound

        bioactivities_data = [
            (5.0, "CCO"),
            (5.5, "CCCO"),
            (6.0, "CCCCO"),
            (7.5, "c1ccccc1"),
            (8.0, "c1ccncc1"),
            (4.5, "CC(=O)O"),
        ]

        compounds = []
        for pic50_val, smi in bioactivities_data:
            bio = Bioactivity(
                target_name="TestTarget",
                value=pic50_val,
                unit="nM",
                activity_type="IC50",
                pic50=pic50_val,
                value_nm=pic50_val,
            )
            compound = Compound(smiles=smi, bioactivities=[bio])
            compounds.append(compound)

        return CompoundCollection(compounds=compounds)

    def test_returns_dict_with_coords(self):
        """plot_sar_landscape returns dict with 'coords' key."""
        collection = self._make_collection_with_bioactivities()

        with (
            patch("chemfuse.analyze.landscape.reduce_dimensions", side_effect=_mock_reduce_dimensions),
            patch("chemfuse.analyze.landscape._detect_cliff_pairs", return_value=[]),
            patch("chemfuse.analyze.landscape._build_plotly_figure", return_value=None),
        ):
            result = collection.plot_sar_landscape(activity_col="pic50", method="pca")

        assert "coords" in result
        assert isinstance(result["coords"], np.ndarray)
        assert result["coords"].shape[1] == 2

    def test_no_valid_activities_raises_value_error(self):
        """Raises ValueError when no compounds have valid activity data."""
        from chemfuse.models.collection import CompoundCollection
        from chemfuse.models.compound import Compound

        collection = CompoundCollection(
            compounds=[
                Compound(smiles="CCO"),
                Compound(smiles="CCCO"),
            ]
        )
        with pytest.raises(ValueError, match="No compounds with valid"):
            collection.plot_sar_landscape(activity_col="pic50", method="pca")

    def test_show_cliffs_false_returns_empty_cliff_pairs(self):
        """show_cliffs=False results in empty cliff_pairs."""
        collection = self._make_collection_with_bioactivities()

        with (
            patch("chemfuse.analyze.landscape.reduce_dimensions", side_effect=_mock_reduce_dimensions),
            patch("chemfuse.analyze.landscape._build_plotly_figure", return_value=None),
        ):
            result = collection.plot_sar_landscape(
                activity_col="pic50",
                method="pca",
                show_cliffs=False,
            )

        # With show_cliffs=False the threshold > 1.0 so no pairs found
        assert result["cliff_pairs"] == []

    def test_coords_length_matches_valid_compounds(self):
        """Number of coordinate rows matches number of compounds with activities."""
        collection = self._make_collection_with_bioactivities()

        with (
            patch("chemfuse.analyze.landscape.reduce_dimensions", side_effect=_mock_reduce_dimensions),
            patch("chemfuse.analyze.landscape._detect_cliff_pairs", return_value=[]),
            patch("chemfuse.analyze.landscape._build_plotly_figure", return_value=None),
        ):
            result = collection.plot_sar_landscape(activity_col="pic50", method="pca")

        assert result["coords"].shape[0] == len(collection.compounds)
