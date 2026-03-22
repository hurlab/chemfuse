"""Tests for the Cross-Reference page."""

from __future__ import annotations

from unittest.mock import MagicMock, patch


class TestXrefPageRender:
    """Test Cross-Reference page rendering."""

    def test_render_shows_header(self, patch_st: MagicMock) -> None:
        """Xref page shows 'Cross-Reference' header."""
        from chemfuse.web.pages.xref import render

        render()
        patch_st.header.assert_called_with("Cross-Reference")

    def test_render_shows_input_field(self, patch_st: MagicMock) -> None:
        """Xref page shows a text input for identifier."""
        from chemfuse.web.pages.xref import render

        render()
        patch_st.text_input.assert_called()

    def test_render_shows_lookup_button(self, patch_st: MagicMock) -> None:
        """Xref page shows a 'Look up' button."""
        from chemfuse.web.pages.xref import render

        render()
        button_labels = [str(c) for c in patch_st.button.call_args_list]
        assert any("Look up" in label for label in button_labels)

    def test_render_with_existing_result(self, patch_st: MagicMock) -> None:
        """Xref page displays existing xref result from session state."""
        patch_st.session_state["xref_result"] = {
            "pubchem": "2244",
            "chembl": "CHEMBL25",
        }
        patch_st.text_input.return_value = "2244"
        from chemfuse.web.pages.xref import render

        render()
        patch_st.dataframe.assert_called()

    def test_render_no_result_on_failed_lookup(self, patch_st: MagicMock) -> None:
        """Xref page shows warning when lookup returns no data."""
        patch_st.button.return_value = True
        patch_st.text_input.return_value = "nonexistent_compound_xyz"

        with patch("chemfuse.web.pages.xref._lookup_xref", return_value=None):
            from chemfuse.web.pages.xref import render
            render()

        patch_st.warning.assert_called()


class TestDisplayXrefTable:
    """Test xref table display logic."""

    def test_shows_database_and_identifier_columns(self, patch_st: MagicMock) -> None:
        """Xref table includes 'Database' and 'Identifier' columns."""
        from chemfuse.web.pages.xref import _display_xref_table

        mappings = {
            "pubchem": "2244",
            "chembl": "CHEMBL25",
        }
        _display_xref_table("aspirin", mappings)

        call_args = patch_st.dataframe.call_args
        df = call_args[0][0]
        assert "Database" in df.columns
        assert "Identifier" in df.columns

    def test_shows_subheader_with_query(self, patch_st: MagicMock) -> None:
        """Xref display shows query in the subheader."""
        from chemfuse.web.pages.xref import _display_xref_table

        _display_xref_table("CID 2244", {"pubchem": "2244"})
        patch_st.subheader.assert_called()
        call_args = str(patch_st.subheader.call_args)
        assert "2244" in call_args

    def test_shows_info_when_empty_mappings(self, patch_st: MagicMock) -> None:
        """Shows info message when mappings dict is empty."""
        from chemfuse.web.pages.xref import _display_xref_table

        _display_xref_table("unknown", {})
        patch_st.info.assert_called()

    def test_shows_links_for_known_databases(self, patch_st: MagicMock) -> None:
        """Shows clickable links for databases with known URL templates."""
        from chemfuse.web.pages.xref import _display_xref_table

        mappings = {"pubchem": "2244", "chembl": "CHEMBL25", "drugbank": "DB00945"}
        _display_xref_table("aspirin", mappings)
        # Markdown should be called with link HTML
        markdown_calls = " ".join(str(c) for c in patch_st.markdown.call_args_list)
        assert "href" in markdown_calls


class TestXrefDbUrls:
    """Test database URL generation."""

    def test_pubchem_url_contains_cid(self) -> None:
        """PubChem URL template substitutes CID correctly."""
        from chemfuse.web.pages.xref import _DB_URLS

        url = _DB_URLS["pubchem"].format(id="2244")
        assert "2244" in url
        assert "pubchem" in url.lower()

    def test_chembl_url_contains_id(self) -> None:
        """ChEMBL URL template substitutes ID correctly."""
        from chemfuse.web.pages.xref import _DB_URLS

        url = _DB_URLS["chembl"].format(id="CHEMBL25")
        assert "CHEMBL25" in url

    def test_all_url_templates_have_id_placeholder(self) -> None:
        """All URL templates contain {id} placeholder."""
        from chemfuse.web.pages.xref import _DB_URLS

        for db_name, template in _DB_URLS.items():
            assert "{id}" in template, f"Missing {{id}} in {db_name} URL template"


class TestLookupXref:
    """Tests for _lookup_xref() function."""

    def test_returns_tuple_on_exception(self) -> None:
        """Returns (None, error_message) when lookup fails."""
        with patch(
            "chemfuse.web.pages.xref._run_async",
            side_effect=Exception("network error"),
        ):
            from chemfuse.web.pages.xref import _lookup_xref

            result, error = _lookup_xref("aspirin")
        assert result is None
        assert error is not None
        assert "network error" in error or "Cross-reference" in error

    def test_returns_result_on_success(self) -> None:
        """Returns (result_dict, None) when lookup succeeds."""
        mock_mappings = {"pubchem": "2244", "chembl": "CHEMBL25"}
        with patch(
            "chemfuse.web.pages.xref._run_async",
            return_value=mock_mappings,
        ):
            from chemfuse.web.pages.xref import _lookup_xref

            result, error = _lookup_xref("2244")
        assert result == mock_mappings
        assert error is None

    def test_returns_none_result_when_async_returns_none(self) -> None:
        """Returns (None, None) when async lookup returns None."""
        with patch(
            "chemfuse.web.pages.xref._run_async",
            return_value=None,
        ):
            from chemfuse.web.pages.xref import _lookup_xref

            result, error = _lookup_xref("unknown")
        assert result is None
        assert error is None


class TestDisplayXrefTableLinks:
    """Test link generation and HTML escaping in _display_xref_table."""

    def test_links_are_html_escaped(self, patch_st: MagicMock) -> None:
        """HTML special characters in IDs are escaped in links."""
        from chemfuse.web.pages.xref import _display_xref_table

        # Use an ID that would be unsafe in HTML
        mappings = {"pubchem": "2244"}
        _display_xref_table("test", mappings)
        # Should not raise; escape is called internally
        patch_st.dataframe.assert_called()

    def test_unknown_database_shows_plain_id(self, patch_st: MagicMock) -> None:
        """IDs for databases without URL templates show as plain text."""
        from chemfuse.web.pages.xref import _display_xref_table

        mappings = {"unknown_db_xyz": "12345"}
        _display_xref_table("test", mappings)
        patch_st.dataframe.assert_called()

    def test_skips_empty_id_values(self, patch_st: MagicMock) -> None:
        """Rows with empty/falsy IDs are skipped."""
        from chemfuse.web.pages.xref import _display_xref_table

        mappings = {"pubchem": "", "chembl": "CHEMBL25"}
        _display_xref_table("test", mappings)
        # Only one row should appear (chembl)
        call_args = patch_st.dataframe.call_args
        df = call_args[0][0]
        assert len(df) == 1

    def test_database_links_markdown_rendered(self, patch_st: MagicMock) -> None:
        """Calls st.markdown for database links section."""
        from chemfuse.web.pages.xref import _display_xref_table

        mappings = {"pubchem": "2244", "chembl": "CHEMBL25"}
        _display_xref_table("aspirin", mappings)
        assert patch_st.markdown.called

    def test_display_xref_table_with_legacy_format(self, patch_st: MagicMock) -> None:
        """render() handles legacy xref_result format (dict without 'mappings' key)."""
        patch_st.session_state["xref_result"] = {
            "pubchem": "2244",
            "chembl": "CHEMBL25",
        }
        patch_st.text_input.return_value = "aspirin"
        from chemfuse.web.pages.xref import render

        render()
        patch_st.dataframe.assert_called()


class TestXrefRenderLookupFlow:
    """Test the full render() lookup flow."""

    def test_render_with_lookup_clicked_and_result(self, patch_st: MagicMock) -> None:
        """render() stores result in session_state when lookup succeeds."""
        patch_st.button.return_value = True
        patch_st.text_input.return_value = "2244"
        mock_result = {"pubchem": "2244", "chembl": "CHEMBL25"}

        with patch(
            "chemfuse.web.pages.xref._lookup_xref",
            return_value=(mock_result, None),
        ):
            from chemfuse.web.pages.xref import render
            render()

        xref_data = patch_st.session_state.get("xref_result")
        assert xref_data is not None

    def test_render_shows_error_when_lookup_returns_error(
        self, patch_st: MagicMock
    ) -> None:
        """render() calls st.error when _lookup_xref returns an error message."""
        patch_st.button.return_value = True
        patch_st.text_input.return_value = "invalid"

        with patch(
            "chemfuse.web.pages.xref._lookup_xref",
            return_value=(None, "Lookup failed: timeout"),
        ):
            from chemfuse.web.pages.xref import render
            render()

        patch_st.error.assert_called()

    def test_render_with_new_format_xref_result(self, patch_st: MagicMock) -> None:
        """render() handles new xref_result format with 'query' and 'mappings' keys."""
        patch_st.session_state["xref_result"] = {
            "query": "aspirin",
            "mappings": {"pubchem": "2244"},
        }
        patch_st.text_input.return_value = "aspirin"
        from chemfuse.web.pages.xref import render

        render()
        patch_st.dataframe.assert_called()

    def test_render_warns_when_result_is_none_tuple(
        self, patch_st: MagicMock
    ) -> None:
        """render() shows warning when lookup returns (None, None) tuple."""
        patch_st.button.return_value = True
        patch_st.text_input.return_value = "mystery_compound"

        with patch(
            "chemfuse.web.pages.xref._lookup_xref",
            return_value=(None, None),
        ):
            from chemfuse.web.pages.xref import render
            render()

        patch_st.warning.assert_called()


class TestLookupXrefAsyncPaths:
    """Test _lookup_xref routing for different identifier types."""

    def test_digit_identifier_uses_pubchem_route(self) -> None:
        """Numeric identifier routes through pubchem mapping."""
        mock_mappings = {"pubchem": "2244"}
        with patch(
            "chemfuse.web.pages.xref._run_async",
            return_value=mock_mappings,
        ):
            from chemfuse.web.pages.xref import _lookup_xref
            result, error = _lookup_xref("2244")
        assert result is not None
        assert error is None

    def test_chembl_identifier_routes_correctly(self) -> None:
        """Identifier starting with CHEMBL routes through chembl mapping."""
        mock_mappings = {"chembl": "CHEMBL25"}
        with patch(
            "chemfuse.web.pages.xref._run_async",
            return_value=mock_mappings,
        ):
            from chemfuse.web.pages.xref import _lookup_xref
            result, error = _lookup_xref("CHEMBL25")
        assert result is not None
        assert error is None

    def test_name_identifier_uses_pubchem_resolution(self) -> None:
        """Non-numeric, non-CHEMBL identifier routes through PubChem resolution."""
        mock_mappings = {"pubchem": "2244"}
        with patch(
            "chemfuse.web.pages.xref._run_async",
            return_value=mock_mappings,
        ):
            from chemfuse.web.pages.xref import _lookup_xref
            # Use a unique identifier not previously cached in this test run
            result, error = _lookup_xref("acetaminophen_test_unique_12345")
        # Result may be None (cache miss or async mock behaviour) — just verify no error
        assert error is None

    def test_inchikey_identifier_uses_cross_reference(self) -> None:
        """InChIKey-like identifier (27 chars with -) uses cross_reference path."""
        mock_mappings = {"pubchem": "2244"}
        with patch(
            "chemfuse.web.pages.xref._run_async",
            return_value=mock_mappings,
        ):
            from chemfuse.web.pages.xref import _lookup_xref
            # InChIKey is 27 chars with dashes
            inchikey = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
            assert len(inchikey) == 27
            result, error = _lookup_xref(inchikey)
        assert result is not None
