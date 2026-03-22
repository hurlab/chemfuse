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

    def test_log_transform_produces_different_scores_than_linear(self):
        """log_transform=True yields different cliff_score than linear for nM activities."""
        with patch("chemfuse.analyze.sar.tanimoto_matrix") as mock_mat:
            mock_mat.return_value = np.array([
                [1.0, 0.9],
                [0.9, 1.0],
            ])

            # 1 nM vs 10000 nM: linear diff = 9999, pIC50 diff ~ 4.0
            compounds = [
                {"smiles": ASPIRIN, "activity": 1.0},
                {"smiles": CAFFEINE, "activity": 10000.0},
            ]

            from chemfuse.analyze.sar import detect_activity_cliffs

            cliffs_linear = detect_activity_cliffs(
                compounds, activity_col="activity", sim_threshold=0.8, log_transform=False
            )
            cliffs_log = detect_activity_cliffs(
                compounds, activity_col="activity", sim_threshold=0.8, log_transform=True
            )

        assert len(cliffs_linear) == 1
        assert len(cliffs_log) == 1
        # Linear score is dominated by the raw nM scale; log score is much smaller
        assert cliffs_linear[0]["cliff_score"] != pytest.approx(cliffs_log[0]["cliff_score"])
        assert cliffs_linear[0]["cliff_score"] > cliffs_log[0]["cliff_score"]

    def test_log_transform_activity_diff_reflects_pic50(self):
        """With log_transform, activity_diff is in pIC50 units, not raw nM."""
        import math

        with patch("chemfuse.analyze.sar.tanimoto_matrix") as mock_mat:
            mock_mat.return_value = np.array([
                [1.0, 0.9],
                [0.9, 1.0],
            ])

            # 10 nM vs 100000 nM -> pIC50 diff = log10(100000/10) = 4.0
            compounds = [
                {"smiles": ASPIRIN, "activity": 10.0},
                {"smiles": CAFFEINE, "activity": 100000.0},
            ]

            from chemfuse.analyze.sar import detect_activity_cliffs
            cliffs = detect_activity_cliffs(
                compounds, activity_col="activity", sim_threshold=0.8, log_transform=True
            )

        assert len(cliffs) == 1
        expected_diff = abs(
            -math.log10(10.0 * 1e-9) - (-math.log10(100000.0 * 1e-9))
        )
        assert cliffs[0]["activity_diff"] == pytest.approx(expected_diff, abs=1e-3)

    def test_log_transform_raw_activities_preserved_in_output(self):
        """With log_transform, activity_1 and activity_2 still hold original nM values."""
        with patch("chemfuse.analyze.sar.tanimoto_matrix") as mock_mat:
            mock_mat.return_value = np.array([
                [1.0, 0.9],
                [0.9, 1.0],
            ])

            compounds = [
                {"smiles": ASPIRIN, "activity": 50.0},
                {"smiles": CAFFEINE, "activity": 5000.0},
            ]

            from chemfuse.analyze.sar import detect_activity_cliffs
            cliffs = detect_activity_cliffs(
                compounds, activity_col="activity", sim_threshold=0.8, log_transform=True
            )

        assert len(cliffs) == 1
        assert cliffs[0]["activity_1"] == pytest.approx(50.0)
        assert cliffs[0]["activity_2"] == pytest.approx(5000.0)

    def test_sali_score_high_sim_high_diff_gives_high_score(self):
        """SALI: high similarity + large activity difference produces a high score."""
        with patch("chemfuse.analyze.sar.tanimoto_matrix") as mock_mat:
            # Very high similarity (0.99) amplifies the SALI denominator (1-0.99 = 0.01)
            mock_mat.return_value = np.array([
                [1.0, 0.99],
                [0.99, 1.0],
            ])

            compounds = [
                {"smiles": ASPIRIN, "activity": 1.0},
                {"smiles": CAFFEINE, "activity": 1001.0},
            ]

            from chemfuse.analyze.sar import detect_activity_cliffs
            cliffs = detect_activity_cliffs(
                compounds, activity_col="activity", sim_threshold=0.8, score_method="sali"
            )

        assert len(cliffs) == 1
        # SALI = |1001 - 1| / (1 - 0.99) = 1000 / 0.01 = 100000
        expected = round(1000.0 / (1.0 - 0.99), 4)
        assert cliffs[0]["cliff_score"] == pytest.approx(expected, rel=1e-3)

    def test_sali_score_greater_than_product_for_high_similarity(self):
        """For high-similarity pairs, SALI score exceeds the product score."""
        with patch("chemfuse.analyze.sar.tanimoto_matrix") as mock_mat:
            mock_mat.return_value = np.array([
                [1.0, 0.95],
                [0.95, 1.0],
            ])

            compounds = [
                {"smiles": ASPIRIN, "activity": 10.0},
                {"smiles": CAFFEINE, "activity": 110.0},
            ]

            from chemfuse.analyze.sar import detect_activity_cliffs

            cliffs_product = detect_activity_cliffs(
                compounds, activity_col="activity", sim_threshold=0.8, score_method="product"
            )
            cliffs_sali = detect_activity_cliffs(
                compounds, activity_col="activity", sim_threshold=0.8, score_method="sali"
            )

        assert len(cliffs_product) == 1
        assert len(cliffs_sali) == 1
        # SALI = 100 / 0.05 = 2000; product = 0.95 * 100 = 95
        assert cliffs_sali[0]["cliff_score"] > cliffs_product[0]["cliff_score"]

    def test_invalid_score_method_raises(self):
        """Unknown score_method raises ValueError."""
        from chemfuse.analyze.sar import detect_activity_cliffs

        compounds = [
            {"smiles": ASPIRIN, "activity": 10.0},
            {"smiles": CAFFEINE, "activity": 100.0},
        ]
        with pytest.raises(ValueError, match="score_method"):
            detect_activity_cliffs(compounds, score_method="unknown")


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
