"""Tests for the main Streamlit application entry point (web/app.py)."""

from __future__ import annotations

from unittest.mock import MagicMock, patch


class TestInitSessionState:
    """Tests for _init_session_state()."""

    def test_sets_current_page_default(self, patch_st: MagicMock) -> None:
        """Sets 'current_page' to 'Search' when not already set."""
        patch_st.session_state.clear()
        from chemfuse.web.app import _init_session_state

        _init_session_state()
        assert patch_st.session_state["current_page"] == "Search"

    def test_sets_search_results_default(self, patch_st: MagicMock) -> None:
        """Sets 'search_results' to empty list when not already set."""
        patch_st.session_state.clear()
        from chemfuse.web.app import _init_session_state

        _init_session_state()
        assert patch_st.session_state["search_results"] == []

    def test_sets_selected_compound_default(self, patch_st: MagicMock) -> None:
        """Sets 'selected_compound' to None when not already set."""
        patch_st.session_state.clear()
        from chemfuse.web.app import _init_session_state

        _init_session_state()
        assert patch_st.session_state["selected_compound"] is None

    def test_sets_session_compounds_default(self, patch_st: MagicMock) -> None:
        """Sets 'session_compounds' to empty list when not already set."""
        patch_st.session_state.clear()
        from chemfuse.web.app import _init_session_state

        _init_session_state()
        assert patch_st.session_state["session_compounds"] == []

    def test_sets_batch_results_default(self, patch_st: MagicMock) -> None:
        """Sets 'batch_results' to None when not already set."""
        patch_st.session_state.clear()
        from chemfuse.web.app import _init_session_state

        _init_session_state()
        assert patch_st.session_state["batch_results"] is None

    def test_sets_xref_result_default(self, patch_st: MagicMock) -> None:
        """Sets 'xref_result' to None when not already set."""
        patch_st.session_state.clear()
        from chemfuse.web.app import _init_session_state

        _init_session_state()
        assert patch_st.session_state["xref_result"] is None

    def test_does_not_overwrite_existing_values(self, patch_st: MagicMock) -> None:
        """Does not overwrite existing session state values."""
        patch_st.session_state["current_page"] = "Profile"
        patch_st.session_state["search_results"] = [{"name": "aspirin"}]
        from chemfuse.web.app import _init_session_state

        _init_session_state()
        assert patch_st.session_state["current_page"] == "Profile"
        assert patch_st.session_state["search_results"] == [{"name": "aspirin"}]

    def test_all_default_keys_present(self, patch_st: MagicMock) -> None:
        """All expected default keys are set."""
        patch_st.session_state.clear()
        from chemfuse.web.app import _init_session_state

        _init_session_state()
        expected_keys = {
            "current_page",
            "search_results",
            "selected_compound",
            "session_compounds",
            "batch_results",
            "xref_result",
        }
        for key in expected_keys:
            assert key in patch_st.session_state, f"Missing key: {key}"


class TestRenderSidebar:
    """Tests for _render_sidebar().

    The sidebar uses ``with st.sidebar:`` context manager, so calls inside
    the block go to the *context manager object* returned by _make_ctx(),
    not to patch_st.sidebar directly. We verify via st.radio and st.title
    since the patch replaces the entire streamlit module.
    """

    def test_renders_without_error(self, patch_st: MagicMock) -> None:
        """_render_sidebar() runs without raising exceptions."""
        patch_st.session_state["current_page"] = "Search"
        patch_st.radio.return_value = "Search"
        from chemfuse.web.app import _render_sidebar

        # Should not raise
        _render_sidebar()

    def test_calls_radio_with_navigation_label(self, patch_st: MagicMock) -> None:
        """st.radio is called with 'Navigation' label."""
        patch_st.session_state["current_page"] = "Search"
        patch_st.radio.return_value = "Search"
        from chemfuse.web.app import _render_sidebar

        _render_sidebar()
        patch_st.radio.assert_called()
        call_args = patch_st.radio.call_args
        assert call_args[0][0] == "Navigation"

    def test_returns_selected_page(self, patch_st: MagicMock) -> None:
        """Returns the selected page name from st.radio."""
        patch_st.session_state["current_page"] = "Search"
        patch_st.radio.return_value = "Search"
        from chemfuse.web.app import _render_sidebar

        result = _render_sidebar()
        assert result == "Search"

    def test_returns_profile_when_selected(self, patch_st: MagicMock) -> None:
        """Returns 'Profile' when radio returns 'Profile'."""
        patch_st.session_state["current_page"] = "Profile"
        patch_st.radio.return_value = "Profile"
        from chemfuse.web.app import _render_sidebar

        result = _render_sidebar()
        assert result == "Profile"

    def test_updates_session_state_current_page(self, patch_st: MagicMock) -> None:
        """Updates session_state['current_page'] with selected page."""
        patch_st.session_state["current_page"] = "Search"
        patch_st.radio.return_value = "Profile"
        from chemfuse.web.app import _render_sidebar

        _render_sidebar()
        assert patch_st.session_state["current_page"] == "Profile"

    def test_renders_caption_for_selected_page(self, patch_st: MagicMock) -> None:
        """Calls st.caption for the selected page description."""
        patch_st.session_state["current_page"] = "Search"
        patch_st.radio.return_value = "Search"
        from chemfuse.web.app import _render_sidebar

        _render_sidebar()
        assert patch_st.caption.called


class TestRenderPage:
    """Tests for _render_page()."""

    def test_render_search_page(self, patch_st: MagicMock) -> None:
        """Dispatches to search.render() for 'Search' page."""
        with patch("chemfuse.web.pages.search.render") as mock_render:
            from chemfuse.web.app import _render_page
            _render_page("Search")
            mock_render.assert_called_once()

    def test_render_profile_page(self, patch_st: MagicMock) -> None:
        """Dispatches to profile.render() for 'Profile' page."""
        with patch("chemfuse.web.pages.profile.render") as mock_render:
            from chemfuse.web.app import _render_page
            _render_page("Profile")
            mock_render.assert_called_once()

    def test_render_batch_screen_page(self, patch_st: MagicMock) -> None:
        """Dispatches to screen.render() for 'Batch Screen' page."""
        with patch("chemfuse.web.pages.screen.render") as mock_render:
            from chemfuse.web.app import _render_page
            _render_page("Batch Screen")
            mock_render.assert_called_once()

    def test_render_chemical_space_page(self, patch_st: MagicMock) -> None:
        """Dispatches to chemspace.render() for 'Chemical Space' page."""
        with patch("chemfuse.web.pages.chemspace.render") as mock_render:
            from chemfuse.web.app import _render_page
            _render_page("Chemical Space")
            mock_render.assert_called_once()

    def test_render_compare_page(self, patch_st: MagicMock) -> None:
        """Dispatches to compare.render() for 'Compare' page."""
        with patch("chemfuse.web.pages.compare.render") as mock_render:
            from chemfuse.web.app import _render_page
            _render_page("Compare")
            mock_render.assert_called_once()

    def test_render_xref_page(self, patch_st: MagicMock) -> None:
        """Dispatches to xref.render() for 'Cross-Reference' page."""
        with patch("chemfuse.web.pages.xref.render") as mock_render:
            from chemfuse.web.app import _render_page
            _render_page("Cross-Reference")
            mock_render.assert_called_once()

    def test_render_unknown_page_shows_error(self, patch_st: MagicMock) -> None:
        """Shows st.error for unknown page name."""
        from chemfuse.web.app import _render_page

        _render_page("NonexistentPage")
        patch_st.error.assert_called_once()
        error_msg = str(patch_st.error.call_args)
        assert "Unknown page" in error_msg or "NonexistentPage" in error_msg

    def test_render_page_handles_import_exception(self, patch_st: MagicMock) -> None:
        """Shows st.error when page module raises exception."""
        with patch("chemfuse.web.pages.search.render", side_effect=RuntimeError("boom")):
            from chemfuse.web.app import _render_page
            _render_page("Search")
        patch_st.error.assert_called_once()
        error_msg = str(patch_st.error.call_args)
        assert "Search" in error_msg or "boom" in error_msg


class TestMain:
    """Tests for main()."""

    def test_main_calls_set_page_config(self, patch_st: MagicMock) -> None:
        """main() calls st.set_page_config."""
        patch_st.session_state["current_page"] = "Search"
        patch_st.radio.return_value = "Search"
        with patch("chemfuse.web.pages.search.render"):
            from chemfuse.web.app import main
            main()
        patch_st.set_page_config.assert_called_once()

    def test_main_calls_markdown_for_css(self, patch_st: MagicMock) -> None:
        """main() injects CSS via st.markdown."""
        patch_st.session_state["current_page"] = "Search"
        patch_st.radio.return_value = "Search"
        with patch("chemfuse.web.pages.search.render"):
            from chemfuse.web.app import main
            main()
        assert patch_st.markdown.called

    def test_main_page_title_is_chemfuse(self, patch_st: MagicMock) -> None:
        """main() sets page_title to 'ChemFuse'."""
        patch_st.session_state["current_page"] = "Search"
        patch_st.radio.return_value = "Search"
        with patch("chemfuse.web.pages.search.render"):
            from chemfuse.web.app import main
            main()
        call_kwargs = patch_st.set_page_config.call_args[1]
        assert call_kwargs.get("page_title") == "ChemFuse"

    def test_main_layout_is_wide(self, patch_st: MagicMock) -> None:
        """main() configures 'wide' layout."""
        patch_st.session_state["current_page"] = "Search"
        patch_st.radio.return_value = "Search"
        with patch("chemfuse.web.pages.search.render"):
            from chemfuse.web.app import main
            main()
        call_kwargs = patch_st.set_page_config.call_args[1]
        assert call_kwargs.get("layout") == "wide"
