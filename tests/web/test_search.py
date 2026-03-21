"""Tests for the Search page."""

from __future__ import annotations

from unittest.mock import MagicMock, patch


class TestSearchPageRender:
    """Test Search page rendering logic."""

    def test_render_shows_header(self, patch_st: MagicMock) -> None:
        """Search page calls st.header with 'Compound Search'."""
        from chemfuse.web.pages.search import render

        render()
        patch_st.header.assert_called_with("Compound Search")

    def test_render_shows_query_input(self, patch_st: MagicMock) -> None:
        """Search page creates a text_input widget."""
        from chemfuse.web.pages.search import render

        render()
        patch_st.text_input.assert_called()

    def test_render_shows_source_checkboxes(self, patch_st: MagicMock) -> None:
        """Search page creates checkboxes for PubChem and ChEMBL."""
        from chemfuse.web.pages.search import render

        render()
        # Two checkboxes for sources
        calls = [str(c) for c in patch_st.checkbox.call_args_list]
        sources = " ".join(calls)
        assert "PubChem" in sources or "pubchem" in sources.lower()

    def test_render_shows_search_type_selector(self, patch_st: MagicMock) -> None:
        """Search page creates a selectbox for search type."""
        from chemfuse.web.pages.search import render

        render()
        patch_st.selectbox.assert_called()

    def test_render_with_existing_results(
        self, patch_st: MagicMock, compound_list: list[dict]
    ) -> None:
        """Search page displays existing session results without re-searching."""
        patch_st.session_state["search_results"] = compound_list
        from chemfuse.web.pages.search import render

        render()
        # Should call subheader or dataframe for results
        called = patch_st.subheader.called or patch_st.dataframe.called
        assert called

    def test_no_search_without_query(self, patch_st: MagicMock) -> None:
        """Search does not fire when query is empty string."""
        patch_st.button.return_value = True
        patch_st.text_input.return_value = ""
        from chemfuse.web.pages.search import render

        render()
        # Should not call spinner (no search executed)
        patch_st.spinner.assert_not_called()

    def test_search_warns_with_no_sources(self, patch_st: MagicMock) -> None:
        """Search warns when no source is selected."""
        patch_st.button.return_value = True
        patch_st.text_input.return_value = "aspirin"
        # Both source checkboxes return False
        patch_st.checkbox.return_value = False
        from chemfuse.web.pages.search import render

        render()
        patch_st.warning.assert_called()


class TestSearchResultsDisplay:
    """Test result display logic."""

    def test_display_results_table_mode(
        self, patch_st: MagicMock, compound_list: list[dict]
    ) -> None:
        """Results display as dataframe in Table mode."""
        patch_st.radio.return_value = "Table"
        patch_st.session_state["search_results"] = compound_list
        from chemfuse.web.pages.search import render

        render()
        patch_st.dataframe.assert_called()

    def test_display_results_grid_mode(
        self, patch_st: MagicMock, compound_list: list[dict]
    ) -> None:
        """Results display as molecule grid in Grid mode."""
        patch_st.radio.return_value = "Grid"
        patch_st.session_state["search_results"] = compound_list
        from chemfuse.web.pages.search import render

        with patch("chemfuse.web.components.mol_grid.render_mol_grid") as mock_grid:
            render()
            mock_grid.assert_called()


class TestSearchExport:
    """Test CSV export functionality."""

    def test_export_csv_generates_download(
        self, patch_st: MagicMock, compound_list: list[dict]
    ) -> None:
        """Export CSV button triggers st.download_button."""
        patch_st.session_state["search_results"] = compound_list
        patch_st.button.side_effect = lambda label, **kw: label == "Export CSV"
        from chemfuse.web.pages.search import render

        render()
        patch_st.download_button.assert_called()


class TestCompoundToDict:
    """Tests for _compound_to_dict helper."""

    def test_converts_model_dump(self) -> None:
        """Handles objects with model_dump method."""
        from chemfuse.web.pages.search import _compound_to_dict

        obj = MagicMock()
        obj.to_dict = None
        del obj.to_dict
        obj.model_dump.return_value = {"name": "test"}
        result = _compound_to_dict(obj)
        assert result == {"name": "test"}

    def test_converts_to_dict(self) -> None:
        """Handles objects with to_dict method."""
        from chemfuse.web.pages.search import _compound_to_dict

        obj = MagicMock()
        obj.to_dict.return_value = {"name": "aspirin"}
        result = _compound_to_dict(obj)
        assert result == {"name": "aspirin"}

    def test_none_returns_empty(self) -> None:
        """None input returns empty dict."""
        from chemfuse.web.pages.search import _compound_to_dict

        result = _compound_to_dict(None)
        assert result == {}
