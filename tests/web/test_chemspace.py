"""Tests for the Chemical Space page."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pandas as pd


class TestChemSpacePageRender:
    """Test Chemical Space page rendering."""

    def test_render_shows_header(self, patch_st: MagicMock) -> None:
        """Chemical Space page shows proper header."""
        from chemfuse.web.pages.chemspace import render

        render()
        patch_st.header.assert_called_with("Chemical Space Visualization")

    def test_render_shows_data_source_selector(self, patch_st: MagicMock) -> None:
        """Chemical Space page creates a radio selector for data source."""
        from chemfuse.web.pages.chemspace import render

        render()
        patch_st.radio.assert_called()

    def test_render_shows_method_selector(self, patch_st: MagicMock) -> None:
        """Chemical Space page shows dimensionality reduction method selector."""
        # Need some compounds in session to show method options
        patch_st.session_state["session_compounds"] = [
            {"smiles": "CCO", "name": "ethanol"},
            {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "name": "aspirin"},
        ]
        patch_st.radio.return_value = "Session results"
        from chemfuse.web.pages.chemspace import render

        render()
        # selectbox called for method
        patch_st.selectbox.assert_called()

    def test_render_shows_info_when_no_session_data(self, patch_st: MagicMock) -> None:
        """Shows info message when session has no compounds."""
        patch_st.session_state["session_compounds"] = []
        patch_st.radio.return_value = "Session results"
        from chemfuse.web.pages.chemspace import render

        render()
        patch_st.info.assert_called()

    def test_render_from_upload(self, patch_st: MagicMock) -> None:
        """Chemical Space page accepts uploaded CSV."""
        patch_st.radio.return_value = "Upload CSV"
        df = pd.DataFrame({"smiles": ["CCO", "CCCO", "CC(=O)O"]})
        mock_file = MagicMock()
        patch_st.file_uploader.return_value = mock_file
        with patch("pandas.read_csv", return_value=df):
            from chemfuse.web.pages.chemspace import render
            render()
        patch_st.success.assert_called()


class TestFindSmilesColumn:
    """Test SMILES column finder in chemspace module."""

    def test_finds_smiles_column(self) -> None:
        from chemfuse.web.pages.chemspace import _find_smiles_column

        df = pd.DataFrame({"smiles": ["CCO"], "mw": [46.0]})
        assert _find_smiles_column(df) == "smiles"

    def test_case_insensitive(self) -> None:
        from chemfuse.web.pages.chemspace import _find_smiles_column

        df = pd.DataFrame({"SMILES": ["CCO"]})
        assert _find_smiles_column(df) == "SMILES"

    def test_returns_none_when_missing(self) -> None:
        from chemfuse.web.pages.chemspace import _find_smiles_column

        df = pd.DataFrame({"name": ["ethanol"]})
        assert _find_smiles_column(df) is None


class TestComputeAndPlot:
    """Test chemical space computation and plotting."""

    def test_too_few_smiles_shows_warning(self, patch_st: MagicMock) -> None:
        """Shows warning when fewer than 3 SMILES provided."""
        df = pd.DataFrame({"smiles": ["CCO", "CCCO"]})
        from chemfuse.web.pages.chemspace import _compute_and_plot

        _compute_and_plot(df, "smiles", "PCA", "morgan", "cluster")
        patch_st.warning.assert_called()

    def test_successful_computation_stores_plot_data(self, patch_st: MagicMock) -> None:
        """Successful computation stores plot data in session_state."""
        df = pd.DataFrame({
            "smiles": [
                "CCO", "CCCO", "CC(=O)O", "c1ccccc1",
                "CC(=O)Oc1ccccc1C(=O)O",
            ],
            "name": ["a", "b", "c", "d", "e"],
        })

        fake_coords = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0], [9.0, 10.0]]
        with patch("chemfuse.analyze.chemspace.reduce_dimensions", return_value=fake_coords):
            from chemfuse.web.pages.chemspace import _compute_and_plot
            _compute_and_plot(df, "smiles", "PCA", "morgan", "cluster")

        assert patch_st.session_state.get("chemspace_plot_data") is not None

    def test_computation_failure_shows_error(self, patch_st: MagicMock) -> None:
        """Shows error when chemspace computation fails."""
        df = pd.DataFrame({
            "smiles": ["CCO", "CCCO", "CC(=O)O"],
        })
        with patch("chemfuse.analyze.chemspace.reduce_dimensions", side_effect=Exception("fail")):
            from chemfuse.web.pages.chemspace import _compute_and_plot
            _compute_and_plot(df, "smiles", "PCA", "morgan", "cluster")

        patch_st.error.assert_called()


class TestRenderPlot:
    """Test Plotly chart rendering."""

    def test_renders_plotly_chart(self, patch_st: MagicMock) -> None:
        """Renders a plotly chart for valid plot data."""

        df = pd.DataFrame({
            "x": [1.0, 2.0, 3.0],
            "y": [4.0, 5.0, 6.0],
            "name": ["a", "b", "c"],
            "cluster": [0, 1, 0],
        })
        patch_st.session_state["chemspace_method"] = "PCA"
        from chemfuse.web.pages.chemspace import _render_plot

        _render_plot(df, "cluster")
        patch_st.plotly_chart.assert_called()

    def test_skips_chart_without_xy_columns(self, patch_st: MagicMock) -> None:
        """Skips rendering when x/y columns are missing."""
        df = pd.DataFrame({"name": ["a", "b"]})
        from chemfuse.web.pages.chemspace import _render_plot

        _render_plot(df, "cluster")
        patch_st.plotly_chart.assert_not_called()
