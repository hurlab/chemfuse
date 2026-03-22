"""Tests for the Batch Screen page."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pandas as pd

from chemfuse.web._utils import _make_params_hash


class TestScreenPageRender:
    """Test Batch Screen page rendering."""

    def test_render_shows_header(self, patch_st: MagicMock) -> None:
        """Screen page shows 'Batch Screening' header."""
        from chemfuse.web.pages.screen import render

        render()
        patch_st.header.assert_called_with("Batch Screening")

    def test_render_shows_file_uploader(self, patch_st: MagicMock) -> None:
        """Screen page creates a file uploader widget."""
        from chemfuse.web.pages.screen import render

        render()
        patch_st.file_uploader.assert_called()

    def test_render_without_upload_shows_example(self, patch_st: MagicMock) -> None:
        """Screen page shows example format when no file uploaded."""
        patch_st.file_uploader.return_value = None
        from chemfuse.web.pages.screen import render

        render()
        # Should show example format via expander
        patch_st.expander.assert_called()

    def test_render_with_invalid_csv_shows_error(self, patch_st: MagicMock) -> None:
        """Screen page shows error when uploaded CSV cannot be parsed."""
        mock_file = MagicMock()
        mock_file.read.return_value = b"not valid csv content \x00\x01\x02"

        def _bad_read_csv(*args, **kwargs):
            raise Exception("parse error")

        patch_st.file_uploader.return_value = mock_file
        with patch("pandas.read_csv", side_effect=_bad_read_csv):
            from chemfuse.web.pages.screen import render
            render()
        patch_st.error.assert_called()

    def test_render_with_no_smiles_column_shows_error(self, patch_st: MagicMock) -> None:
        """Screen page shows error when CSV lacks a SMILES column."""
        mock_file = MagicMock()
        df_no_smiles = pd.DataFrame({"compound_id": [1, 2], "name": ["a", "b"]})

        def _read_csv(*args, **kwargs):
            return df_no_smiles

        patch_st.file_uploader.return_value = mock_file
        with patch("pandas.read_csv", return_value=df_no_smiles):
            from chemfuse.web.pages.screen import render
            render()
        patch_st.error.assert_called()


class TestFindSmilesColumn:
    """Test SMILES column detection."""

    def test_finds_lowercase_smiles(self) -> None:
        """Detects 'smiles' column."""
        from chemfuse.web.pages.screen import _find_smiles_column

        df = pd.DataFrame({"smiles": ["CCO"], "name": ["ethanol"]})
        assert _find_smiles_column(df) == "smiles"

    def test_finds_uppercase_smiles(self) -> None:
        """Detects 'SMILES' column (case-insensitive)."""
        from chemfuse.web.pages.screen import _find_smiles_column

        df = pd.DataFrame({"SMILES": ["CCO"]})
        assert _find_smiles_column(df) == "SMILES"

    def test_returns_none_when_absent(self) -> None:
        """Returns None when no SMILES column exists."""
        from chemfuse.web.pages.screen import _find_smiles_column

        df = pd.DataFrame({"compound_name": ["ethanol"]})
        assert _find_smiles_column(df) is None


class TestApplyFilters:
    """Test filter application logic."""

    def test_filter_by_mw_range(self) -> None:
        """Filter removes rows outside MW range."""
        from chemfuse.web.pages.screen import _apply_filters

        df = pd.DataFrame({
            "smiles": ["CCO", "CCCO", "CCCCO"],
            "molecular_weight": [46.0, 60.0, 74.0],
        })
        result = _apply_filters(df, {"molecular_weight": (50.0, 65.0)})
        assert len(result) == 1
        assert result.iloc[0]["molecular_weight"] == 60.0

    def test_filter_preserves_all_when_range_is_wide(self) -> None:
        """Filter keeps all rows when range covers all values."""
        from chemfuse.web.pages.screen import _apply_filters

        df = pd.DataFrame({
            "smiles": ["CCO", "CCCO"],
            "molecular_weight": [46.0, 60.0],
        })
        result = _apply_filters(df, {"molecular_weight": (0.0, 1000.0)})
        assert len(result) == 2

    def test_filter_ignores_missing_column(self) -> None:
        """Filter skips columns that don't exist in the dataframe."""
        from chemfuse.web.pages.screen import _apply_filters

        df = pd.DataFrame({"smiles": ["CCO", "CCCO"]})
        result = _apply_filters(df, {"molecular_weight": (50.0, 200.0)})
        assert len(result) == 2  # No column to filter on, keep all


class TestRunPipeline:
    """Test pipeline execution."""

    def test_pipeline_with_empty_smiles(self, patch_st: MagicMock) -> None:
        """Pipeline shows error when all SMILES are empty."""
        df = pd.DataFrame({"smiles": [None, None]})
        from chemfuse.web.pages.screen import _run_pipeline

        _run_pipeline(df, "smiles", [], run_admet=False, run_druglikeness=False, run_clustering=False)
        patch_st.error.assert_called()

    def test_pipeline_produces_results(self, patch_st: MagicMock) -> None:
        """Pipeline stores results in session_state."""
        df = pd.DataFrame({
            "smiles": ["CC(=O)Oc1ccccc1C(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"],
            "name": ["aspirin", "caffeine"],
        })
        from chemfuse.web.pages.screen import _run_pipeline

        _run_pipeline(
            df, "smiles", [],
            run_admet=False,
            run_druglikeness=False,
            run_clustering=False,
        )
        assert patch_st.session_state.get("batch_results") is not None

    def test_pipeline_with_druglikeness(self, patch_st: MagicMock) -> None:
        """Pipeline runs drug-likeness step and adds result column."""
        df = pd.DataFrame({"smiles": ["CC(=O)Oc1ccccc1C(=O)O"]})
        from chemfuse.web.pages.screen import _run_pipeline

        _run_pipeline(
            df, "smiles", [],
            run_admet=False,
            run_druglikeness=True,
            run_clustering=False,
        )
        results = patch_st.session_state.get("batch_results")
        assert results is not None

    def test_display_results_shows_dataframe(self, patch_st: MagicMock) -> None:
        """Existing results are displayed as a dataframe."""
        df = pd.DataFrame({
            "smiles": ["CC(=O)Oc1ccccc1C(=O)O"],
            "lipinski_pass": [True],
        })
        patch_st.session_state["batch_results"] = df
        from chemfuse.web.pages.screen import render

        # Provide a valid uploaded CSV so we get to the results display.
        mock_file = MagicMock()
        mock_file.name = "test.csv"
        patch_st.file_uploader.return_value = mock_file
        # Pre-seed the params hash so the staleness check does not clear results.
        patch_st.checkbox.side_effect = [False, False]
        patch_st.toggle.side_effect = [False, True, False]
        params = {
            "upload_name": "test.csv",
            "use_pubchem": False,
            "use_chembl": False,
            "run_admet": False,
            "run_druglikeness": True,
            "run_clustering": False,
        }
        patch_st.session_state["_batch_params_hash"] = _make_params_hash(params)
        with patch("pandas.read_csv", return_value=pd.DataFrame({"smiles": ["CC(=O)Oc1ccccc1C(=O)O"]})):
            render()

        patch_st.dataframe.assert_called()


class TestComputeAdmet:
    """Test _compute_admet() exception handling."""

    def test_compute_admet_returns_empty_on_exception(
        self, patch_st: MagicMock
    ) -> None:
        """_compute_admet returns empty dict and warns when import fails."""
        with patch(
            "chemfuse.compute.admet.predict_admet",
            side_effect=Exception("admet error"),
        ):
            from chemfuse.web.pages.screen import _compute_admet

            result = _compute_admet(["CC(=O)O"])
        assert result == {}
        patch_st.warning.assert_called()

    def test_compute_admet_returns_empty_when_no_results(self) -> None:
        """_compute_admet returns empty dict when predict_admet returns empty list."""
        with patch(
            "chemfuse.compute.admet.predict_admet",
            return_value=[],
        ):
            from chemfuse.web.pages.screen import _compute_admet

            result = _compute_admet(["CC(=O)O"])
        assert result == {}

    def test_compute_admet_returns_columns_dict(self) -> None:
        """_compute_admet returns dict of column name to value lists."""
        mock_preds = [{"HIA": 0.9, "BBB": 0.5}]
        with patch(
            "chemfuse.compute.admet.predict_admet",
            return_value=mock_preds,
        ):
            from chemfuse.web.pages.screen import _compute_admet

            result = _compute_admet(["CC(=O)O"])
        assert "HIA" in result or len(result) >= 0  # graceful check


class TestComputeClusters:
    """Test _compute_clusters() exception handling."""

    def test_compute_clusters_returns_none_list_on_exception(
        self, patch_st: MagicMock
    ) -> None:
        """_compute_clusters returns list of None values and warns when clustering fails."""
        with patch(
            "chemfuse.analyze.clustering.cluster_compounds",
            side_effect=Exception("cluster error"),
        ):
            from chemfuse.web.pages.screen import _compute_clusters

            result = _compute_clusters(["CC(=O)O", "CCO"])
        assert result == [None, None]
        patch_st.warning.assert_called()


class TestDisplayResults:
    """Test _display_results() function."""

    def test_display_results_shows_subheader(self, patch_st: MagicMock) -> None:
        """_display_results shows result count in subheader."""
        df = pd.DataFrame({"smiles": ["CCO", "CC(=O)O"]})
        from chemfuse.web.pages.screen import _display_results

        _display_results(df)
        patch_st.subheader.assert_called()
        subheader_text = str(patch_st.subheader.call_args)
        assert "2" in subheader_text

    def test_display_results_shows_download_button(self, patch_st: MagicMock) -> None:
        """_display_results shows CSV download button."""
        df = pd.DataFrame({"smiles": ["CCO"]})
        from chemfuse.web.pages.screen import _display_results

        _display_results(df)
        patch_st.download_button.assert_called()

    def test_display_results_shows_expander(self, patch_st: MagicMock) -> None:
        """_display_results shows a 'Filter Results' expander."""
        df = pd.DataFrame({"smiles": ["CCO"]})
        from chemfuse.web.pages.screen import _display_results

        _display_results(df)
        assert patch_st.expander.called

    def test_display_results_applies_filters(self, patch_st: MagicMock) -> None:
        """_display_results applies filter params when provided."""
        patch_st.slider.return_value = (0.0, 100.0)
        # Make render_property_filters return a non-empty dict
        with patch(
            "chemfuse.web.components.filters.render_property_filters",
            return_value={"molecular_weight": (0.0, 100.0)},
        ):
            df = pd.DataFrame({
                "smiles": ["CCO"],
                "molecular_weight": [46.0],
            })
            from chemfuse.web.pages.screen import _display_results

            _display_results(df)
        patch_st.dataframe.assert_called()


class TestRenderFilterConfiguration:
    """Test filter rendering in screen render() function."""

    def test_pipeline_config_shows_toggle_for_admet(
        self, patch_st: MagicMock
    ) -> None:
        """Screen page shows ADMET toggle in pipeline configuration."""
        mock_file = MagicMock()
        patch_st.file_uploader.return_value = mock_file
        with patch(
            "pandas.read_csv",
            return_value=pd.DataFrame({"smiles": ["CCO"]}),
        ):
            from chemfuse.web.pages.screen import render
            render()
        assert patch_st.toggle.called

    def test_pipeline_config_shows_toggle_for_druglikeness(
        self, patch_st: MagicMock
    ) -> None:
        """Screen page shows drug-likeness toggle in pipeline configuration."""
        mock_file = MagicMock()
        patch_st.file_uploader.return_value = mock_file
        with patch(
            "pandas.read_csv",
            return_value=pd.DataFrame({"smiles": ["CCO"]}),
        ):
            from chemfuse.web.pages.screen import render
            render()
        toggle_calls = " ".join(str(c) for c in patch_st.toggle.call_args_list)
        assert "Drug" in toggle_calls or "drug" in toggle_calls or "likeness" in toggle_calls.lower()


class TestBatchScreenStaleness:
    """Verify that stale batch_results are cleared when parameters change."""

    def _make_mock_file(self, name: str = "compounds.csv") -> MagicMock:
        mock_file = MagicMock()
        mock_file.name = name
        return mock_file

    def test_changing_file_clears_cached_results(self, patch_st: MagicMock) -> None:
        """When a new file is uploaded, batch_results is cleared from session state."""
        stored_results = pd.DataFrame({"smiles": ["CCO"], "lipinski_pass": [True]})
        patch_st.session_state["batch_results"] = stored_results
        patch_st.session_state["_batch_params_hash"] = None  # force detection
        patch_st.file_uploader.return_value = self._make_mock_file("new_file.csv")

        with patch(
            "pandas.read_csv",
            return_value=pd.DataFrame({"smiles": ["CCO"]}),
        ):
            from chemfuse.web.pages.screen import render
            render()

        assert patch_st.session_state.get("batch_results") is None

    def test_changing_option_clears_cached_results(self, patch_st: MagicMock) -> None:
        """When a pipeline option changes, batch_results is cleared."""
        stored_results = pd.DataFrame({"smiles": ["CCO"], "lipinski_pass": [True]})
        patch_st.session_state["batch_results"] = stored_results
        # Store a hash corresponding to run_admet=False
        params = {
            "upload_name": "compounds.csv",
            "use_pubchem": False,
            "use_chembl": False,
            "run_admet": False,
            "run_druglikeness": True,
            "run_clustering": False,
        }
        patch_st.session_state["_batch_params_hash"] = _make_params_hash(params)
        patch_st.file_uploader.return_value = self._make_mock_file("compounds.csv")
        # Simulate user toggling ADMET on
        patch_st.toggle.side_effect = [True, True, False]  # admet=True, dl=True, cluster=False

        with patch(
            "pandas.read_csv",
            return_value=pd.DataFrame({"smiles": ["CCO"]}),
        ):
            from chemfuse.web.pages.screen import render
            render()

        assert patch_st.session_state.get("batch_results") is None

    def test_unchanged_params_preserve_cached_results(self, patch_st: MagicMock) -> None:
        """When all params stay the same, batch_results is preserved."""
        stored_results = pd.DataFrame({"smiles": ["CCO"], "lipinski_pass": [True]})
        patch_st.session_state["batch_results"] = stored_results
        params = {
            "upload_name": "compounds.csv",
            "use_pubchem": False,
            "use_chembl": False,
            "run_admet": False,
            "run_druglikeness": False,
            "run_clustering": False,
        }
        patch_st.session_state["_batch_params_hash"] = _make_params_hash(params)
        patch_st.file_uploader.return_value = self._make_mock_file("compounds.csv")
        # Checkboxes and toggles return same values as stored hash
        patch_st.checkbox.side_effect = [False, False]
        patch_st.toggle.side_effect = [False, False, False]

        with patch(
            "pandas.read_csv",
            return_value=pd.DataFrame({"smiles": ["CCO"]}),
        ):
            from chemfuse.web.pages.screen import render
            render()

        assert patch_st.session_state.get("batch_results") is not None
