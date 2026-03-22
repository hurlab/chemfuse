"""Tests for the Chemical Space page."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pandas as pd

from chemfuse.web._utils import _make_params_hash, _params_changed


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


class TestParamsHashUtility:
    """Tests for the _make_params_hash and _params_changed helpers."""

    def test_same_params_produce_same_hash(self) -> None:
        """Identical parameter dicts produce the same hash."""
        params = {"method": "PCA", "fp_type": "morgan", "color_by": "cluster"}
        assert _make_params_hash(params) == _make_params_hash(params)

    def test_different_params_produce_different_hash(self) -> None:
        """Changing a parameter value changes the hash."""
        params_a = {"method": "PCA", "fp_type": "morgan"}
        params_b = {"method": "UMAP", "fp_type": "morgan"}
        assert _make_params_hash(params_a) != _make_params_hash(params_b)

    def test_params_changed_returns_true_on_first_call(self) -> None:
        """_params_changed returns True when no previous hash is stored."""
        state = {}
        assert _params_changed(state, "_test_hash", {"method": "PCA"}) is True

    def test_params_changed_returns_false_when_params_unchanged(self) -> None:
        """_params_changed returns False when params are identical to last call."""
        state = {}
        params = {"method": "PCA", "fp_type": "morgan"}
        _params_changed(state, "_test_hash", params)  # first call stores hash
        assert _params_changed(state, "_test_hash", params) is False

    def test_params_changed_returns_true_after_param_update(self) -> None:
        """_params_changed returns True after a parameter value changes."""
        state = {}
        _params_changed(state, "_test_hash", {"method": "PCA"})
        assert _params_changed(state, "_test_hash", {"method": "UMAP"}) is True

    def test_params_changed_updates_stored_hash(self) -> None:
        """After a change, the stored hash reflects the new parameter set."""
        state = {}
        _params_changed(state, "_test_hash", {"method": "PCA"})
        _params_changed(state, "_test_hash", {"method": "UMAP"})
        # Now calling again with UMAP should return False (hash was updated)
        assert _params_changed(state, "_test_hash", {"method": "UMAP"}) is False


class TestChemSpaceStaleness:
    """Verify that stale chemspace_plot_data is cleared on parameter change."""

    def _make_df_with_smiles(self) -> pd.DataFrame:
        return pd.DataFrame({
            "smiles": ["CCO", "CCCO", "CC(=O)O", "c1ccccc1"],
            "name": ["a", "b", "c", "d"],
        })

    def test_changing_method_clears_cached_plot(self, patch_st: MagicMock) -> None:
        """When the method changes, chemspace_plot_data is removed from session state."""
        # Seed session state with stale plot data computed with PCA
        patch_st.session_state["chemspace_plot_data"] = pd.DataFrame({"x": [1.0], "y": [2.0]})
        patch_st.session_state["_chemspace_params_hash"] = None  # force detection
        patch_st.radio.return_value = "Session results"
        patch_st.session_state["session_compounds"] = [
            {"smiles": s, "name": n}
            for s, n in zip(
                ["CCO", "CCCO", "CC(=O)O", "c1ccccc1"],
                ["a", "b", "c", "d"],
                strict=True,
            )
        ]
        # Simulate user switching from PCA to UMAP
        patch_st.selectbox.side_effect = ["UMAP", "morgan", "cluster"]

        from chemfuse.web.pages.chemspace import render
        render()

        # Plot data must have been cleared because hash did not match
        assert patch_st.session_state.get("chemspace_plot_data") is None

    def test_unchanged_params_preserve_cached_plot(self, patch_st: MagicMock) -> None:
        """When params are unchanged, chemspace_plot_data is preserved."""
        stored_plot = pd.DataFrame({"x": [1.0], "y": [2.0]})
        patch_st.session_state["chemspace_plot_data"] = stored_plot
        patch_st.radio.return_value = "Session results"
        patch_st.session_state["session_compounds"] = [
            {"smiles": s, "name": n}
            for s, n in zip(
                ["CCO", "CCCO", "CC(=O)O", "c1ccccc1"],
                ["a", "b", "c", "d"],
                strict=True,
            )
        ]
        # Set up a pre-stored hash that matches these params
        from chemfuse.web._utils import _make_params_hash
        params = {
            "source": "Session results",
            "upload_name": None,
            "method": "PCA",
            "fp_type": "morgan",
            "color_by": "cluster",
        }
        patch_st.session_state["_chemspace_params_hash"] = _make_params_hash(params)
        patch_st.selectbox.side_effect = ["PCA", "morgan", "cluster"]

        from chemfuse.web.pages.chemspace import render
        render()

        # Plot data must still be present because params did not change
        assert patch_st.session_state.get("chemspace_plot_data") is not None
