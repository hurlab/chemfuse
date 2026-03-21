"""Tests for the Batch Screen page."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pandas as pd


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

        # Provide a valid uploaded CSV so we get to the results display
        mock_file = MagicMock()
        patch_st.file_uploader.return_value = mock_file
        with patch("pandas.read_csv", return_value=pd.DataFrame({"smiles": ["CC(=O)Oc1ccccc1C(=O)O"]})):
            render()

        patch_st.dataframe.assert_called()
