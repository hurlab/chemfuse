"""Tests for analyze/sar.py (SPEC-CF-005)."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from chemfuse.core.exceptions import OptionalDependencyError

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
CAFFEINE = "Cn1cnc2c1c(=O)n(c(=O)n2C)C"
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"

COMPOUNDS = [
    {"smiles": ASPIRIN, "name": "aspirin", "activity": 10.0},
    {"smiles": CAFFEINE, "name": "caffeine", "activity": 100.0},
    {"smiles": IBUPROFEN, "name": "ibuprofen", "activity": 20.0},
]


class TestDetectActivityCliffs:
    def test_missing_activity_column_raises(self):
        """Raises ValueError when activity column is absent from a compound."""
        from chemfuse.analyze.sar import detect_activity_cliffs

        bad_compounds = [{"smiles": ASPIRIN}, {"smiles": CAFFEINE}]
        with pytest.raises(ValueError, match="activity"):
            detect_activity_cliffs(bad_compounds, activity_col="activity")

    def test_empty_compounds_returns_empty(self):
        """Empty compounds list returns empty cliffs list."""
        from chemfuse.analyze.sar import detect_activity_cliffs
        assert detect_activity_cliffs([]) == []

    def test_returns_cliffs_sorted_by_score(self):
        """detect_activity_cliffs returns cliffs sorted by descending cliff_score."""
        with patch("chemfuse.analyze.sar.tanimoto_matrix") as mock_mat:
            mock_mat.return_value = np.array([
                [1.0, 0.9, 0.2],
                [0.9, 1.0, 0.3],
                [0.2, 0.3, 1.0],
            ])

            from chemfuse.analyze.sar import detect_activity_cliffs
            cliffs = detect_activity_cliffs(COMPOUNDS, activity_col="activity", sim_threshold=0.8)

        assert isinstance(cliffs, list)
        if len(cliffs) > 1:
            scores = [c["cliff_score"] for c in cliffs]
            assert scores == sorted(scores, reverse=True)

    def test_cliff_score_formula(self):
        """cliff_score = similarity * |activity_diff|."""
        with patch("chemfuse.analyze.sar.tanimoto_matrix") as mock_mat:
            mock_mat.return_value = np.array([
                [1.0, 0.9],
                [0.9, 1.0],
            ])

            compounds = [
                {"smiles": ASPIRIN, "activity": 10.0},
                {"smiles": CAFFEINE, "activity": 100.0},
            ]

            from chemfuse.analyze.sar import detect_activity_cliffs
            cliffs = detect_activity_cliffs(compounds, activity_col="activity", sim_threshold=0.8)

        assert len(cliffs) == 1
        cliff = cliffs[0]
        expected_score = round(0.9 * abs(10.0 - 100.0), 4)
        assert cliff["cliff_score"] == pytest.approx(expected_score)

    def test_below_threshold_pairs_excluded(self):
        """Pairs with similarity below threshold are not returned as cliffs."""
        with patch("chemfuse.analyze.sar.tanimoto_matrix") as mock_mat:
            # All pairs below threshold (0.5)
            mock_mat.return_value = np.array([
                [1.0, 0.3],
                [0.3, 1.0],
            ])

            compounds = [
                {"smiles": ASPIRIN, "activity": 10.0},
                {"smiles": CAFFEINE, "activity": 100.0},
            ]

            from chemfuse.analyze.sar import detect_activity_cliffs
            cliffs = detect_activity_cliffs(
                compounds, activity_col="activity", sim_threshold=0.5
            )

        assert len(cliffs) == 0

    def test_high_sim_large_diff_detected(self):
        """Pairs with sim > 0.8 and activity_diff > 100x are detected."""
        with patch("chemfuse.analyze.sar.tanimoto_matrix") as mock_mat:
            mock_mat.return_value = np.array([
                [1.0, 0.85],
                [0.85, 1.0],
            ])

            compounds = [
                {"smiles": ASPIRIN, "activity": 1.0},
                {"smiles": CAFFEINE, "activity": 200.0},
            ]

            from chemfuse.analyze.sar import detect_activity_cliffs
            cliffs = detect_activity_cliffs(
                compounds, activity_col="activity", sim_threshold=0.8
            )

        assert len(cliffs) == 1
        assert cliffs[0]["activity_diff"] == pytest.approx(199.0)

    def test_cliff_dict_has_required_keys(self):
        """Each cliff dict has all required keys."""
        with patch("chemfuse.analyze.sar.tanimoto_matrix") as mock_mat:
            mock_mat.return_value = np.array([
                [1.0, 0.9],
                [0.9, 1.0],
            ])

            compounds = [
                {"smiles": ASPIRIN, "activity": 10.0},
                {"smiles": CAFFEINE, "activity": 100.0},
            ]

            from chemfuse.analyze.sar import detect_activity_cliffs
            cliffs = detect_activity_cliffs(compounds, activity_col="activity", sim_threshold=0.8)

        required_keys = {"smiles_1", "smiles_2", "similarity", "activity_1", "activity_2", "activity_diff", "cliff_score"}
        assert required_keys.issubset(set(cliffs[0].keys()))


def _make_mock_nx():
    """Create a minimal networkx mock with Graph support."""

    mock_nx = MagicMock()

    class MockGraph:
        def __init__(self):
            self._nodes: dict = {}
            self._edges: list = []

        def add_node(self, n: int, **attrs: object) -> None:
            self._nodes[n] = attrs

        def add_edge(self, u: int, v: int, **attrs: object) -> None:
            self._edges.append((u, v, attrs))

        def number_of_nodes(self) -> int:
            return len(self._nodes)

        def number_of_edges(self) -> int:
            return len(self._edges)

        def has_edge(self, u: int, v: int) -> bool:
            return any((a == u and b == v) or (a == v and b == u) for a, b, _ in self._edges)

    mock_nx.Graph = MockGraph
    return mock_nx


class TestBuildSimilarityNetwork:
    def test_requires_networkx(self):
        """Raises OptionalDependencyError when networkx is absent."""
        with patch("chemfuse.analyze.sar._NX_AVAILABLE", False):
            from chemfuse.analyze.sar import build_similarity_network
            with pytest.raises(OptionalDependencyError):
                build_similarity_network(COMPOUNDS)

    def test_returns_graph_with_correct_nodes(self):
        """build_similarity_network returns a graph with n nodes."""
        mock_nx = _make_mock_nx()

        with patch("chemfuse.analyze.sar._NX_AVAILABLE", True), \
             patch("chemfuse.analyze.sar.nx", mock_nx), \
             patch("chemfuse.analyze.sar.tanimoto_matrix") as mock_mat:

            mock_mat.return_value = np.array([
                [1.0, 0.7, 0.3],
                [0.7, 1.0, 0.4],
                [0.3, 0.4, 1.0],
            ])

            from chemfuse.analyze.sar import build_similarity_network
            G = build_similarity_network(COMPOUNDS, threshold=0.5)

        assert G.number_of_nodes() == 3

    def test_graph_edges_respect_threshold(self):
        """Only pairs with similarity >= threshold have edges."""
        mock_nx = _make_mock_nx()

        with patch("chemfuse.analyze.sar._NX_AVAILABLE", True), \
             patch("chemfuse.analyze.sar.nx", mock_nx), \
             patch("chemfuse.analyze.sar.tanimoto_matrix") as mock_mat:

            mock_mat.return_value = np.array([
                [1.0, 0.9, 0.2],
                [0.9, 1.0, 0.3],
                [0.2, 0.3, 1.0],
            ])

            from chemfuse.analyze.sar import build_similarity_network
            G = build_similarity_network(COMPOUNDS, threshold=0.5)

        # Only edge 0-1 should exist (sim=0.9 >= 0.5); others below threshold
        assert G.number_of_edges() == 1
        assert G.has_edge(0, 1)
