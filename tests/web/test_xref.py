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
