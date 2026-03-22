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


class TestExecuteSearch:
    """Tests for _execute_search()."""

    def test_execute_search_calls_spinner(self, patch_st: MagicMock) -> None:
        """_execute_search shows a spinner during search."""
        with patch(
            "chemfuse.web.pages.search._fetch_results",
            return_value=([], []),
        ):
            from chemfuse.web.pages.search import _execute_search

            _execute_search("aspirin", "name", ["pubchem"])
        patch_st.spinner.assert_called()

    def test_execute_search_stores_results(self, patch_st: MagicMock) -> None:
        """_execute_search stores results in session_state."""
        results = [{"name": "aspirin", "smiles": "CC(=O)Oc1ccccc1C(=O)O"}]
        with patch(
            "chemfuse.web.pages.search._fetch_results",
            return_value=(results, []),
        ):
            from chemfuse.web.pages.search import _execute_search

            _execute_search("aspirin", "name", ["pubchem"])
        assert patch_st.session_state.get("search_results") == results

    def test_execute_search_warns_when_no_results(self, patch_st: MagicMock) -> None:
        """_execute_search warns when search returns no results."""
        with patch(
            "chemfuse.web.pages.search._fetch_results",
            return_value=([], []),
        ):
            from chemfuse.web.pages.search import _execute_search

            _execute_search("unknownxyz123", "name", ["pubchem"])
        patch_st.warning.assert_called()

    def test_execute_search_shows_warnings_from_fetch(self, patch_st: MagicMock) -> None:
        """_execute_search shows warnings returned by _fetch_results."""
        warnings = ["Source 'pubchem' error: timeout"]
        with patch(
            "chemfuse.web.pages.search._fetch_results",
            return_value=([], warnings),
        ):
            from chemfuse.web.pages.search import _execute_search

            _execute_search("aspirin", "name", ["pubchem"])
        patch_st.warning.assert_called()

    def test_execute_search_updates_session_compounds(self, patch_st: MagicMock) -> None:
        """_execute_search adds new results to session_compounds."""
        results = [{"name": "aspirin", "smiles": "CC(=O)O"}]
        patch_st.session_state["session_compounds"] = []
        with patch(
            "chemfuse.web.pages.search._fetch_results",
            return_value=(results, []),
        ):
            from chemfuse.web.pages.search import _execute_search

            _execute_search("aspirin", "name", ["pubchem"])
        session_compounds = patch_st.session_state.get("session_compounds", [])
        assert len(session_compounds) > 0

    def test_execute_search_caps_session_compounds(self, patch_st: MagicMock) -> None:
        """_execute_search caps session_compounds at MAX_SESSION_COMPOUNDS."""
        from chemfuse.web.pages.search import MAX_SESSION_COMPOUNDS

        # Pre-fill session with compounds (already at cap)
        existing = [{"name": f"c{i}", "smiles": f"C{i}"} for i in range(MAX_SESSION_COMPOUNDS)]
        patch_st.session_state["session_compounds"] = existing

        new_results = [{"name": "new", "smiles": "CO"}]
        with patch(
            "chemfuse.web.pages.search._fetch_results",
            return_value=(new_results, []),
        ):
            from chemfuse.web.pages.search import _execute_search

            _execute_search("methanol", "name", ["pubchem"])

        session_compounds = patch_st.session_state.get("session_compounds", [])
        assert len(session_compounds) <= MAX_SESSION_COMPOUNDS


class TestFetchResults:
    """Tests for _fetch_results()."""

    def test_returns_empty_on_run_async_error(self) -> None:
        """Returns empty results and a warning when _run_async raises."""
        with patch(
            "chemfuse.web.pages.search._run_async",
            side_effect=Exception("connection refused"),
        ):
            from chemfuse.web.pages.search import _fetch_results

            results, warnings = _fetch_results("aspirin", "name", ("pubchem",))
        assert results == []
        assert len(warnings) > 0

    def test_returns_tuple(self) -> None:
        """Return value is always a tuple of (list, list)."""
        with patch(
            "chemfuse.web.pages.search._run_async",
            side_effect=Exception("err"),
        ):
            from chemfuse.web.pages.search import _fetch_results

            result = _fetch_results("x", "name", ("pubchem",))
        assert isinstance(result, tuple)
        assert len(result) == 2


class TestDisplayResults:
    """Tests for _display_results()."""

    def test_display_results_shows_subheader(
        self, patch_st: MagicMock, compound_list: list[dict]
    ) -> None:
        """_display_results calls st.subheader with result count."""
        patch_st.radio.return_value = "Table"
        from chemfuse.web.pages.search import _display_results

        _display_results(compound_list)
        patch_st.subheader.assert_called()
        subheader_text = str(patch_st.subheader.call_args)
        assert "2" in subheader_text

    def test_display_results_table_mode_calls_dataframe(
        self, patch_st: MagicMock, compound_list: list[dict]
    ) -> None:
        """_display_results calls st.dataframe when in Table view mode."""
        patch_st.radio.return_value = "Table"
        from chemfuse.web.pages.search import _display_results

        _display_results(compound_list)
        patch_st.dataframe.assert_called()

    def test_display_results_grid_mode_calls_mol_grid(
        self, patch_st: MagicMock, compound_list: list[dict]
    ) -> None:
        """_display_results calls render_mol_grid when in Grid view mode."""
        patch_st.radio.return_value = "Grid"
        with patch("chemfuse.web.components.mol_grid.render_mol_grid") as mock_grid:
            from chemfuse.web.pages.search import _display_results

            _display_results(compound_list)
            mock_grid.assert_called()

    def test_display_results_export_csv_button(
        self, patch_st: MagicMock, compound_list: list[dict]
    ) -> None:
        """_display_results renders an 'Export CSV' button."""
        patch_st.radio.return_value = "Table"
        patch_st.button.side_effect = lambda label, **kw: label == "Export CSV"
        from chemfuse.web.pages.search import _display_results

        _display_results(compound_list)
        patch_st.download_button.assert_called()
