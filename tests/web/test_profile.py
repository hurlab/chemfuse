"""Tests for the Compound Profile page."""

from __future__ import annotations

from unittest.mock import MagicMock, patch


class TestProfilePageRender:
    """Test Profile page rendering behavior."""

    def test_render_shows_header(self, patch_st: MagicMock) -> None:
        """Profile page calls st.header."""
        from chemfuse.web.pages.profile import render

        render()
        patch_st.header.assert_called_with("Compound Profile")

    def test_render_with_no_compound_shows_info(self, patch_st: MagicMock) -> None:
        """Profile page shows info message when no compound selected."""
        patch_st.session_state["selected_compound"] = None
        from chemfuse.web.pages.profile import render

        render()
        patch_st.info.assert_called()

    def test_render_with_compound_shows_subheader(
        self, patch_st: MagicMock, aspirin_data: dict
    ) -> None:
        """Profile page shows compound name as subheader."""
        patch_st.session_state["selected_compound"] = aspirin_data
        from chemfuse.web.pages.profile import render

        render()
        patch_st.subheader.assert_called_with("aspirin")

    def test_render_shows_structure_image_for_cid(
        self, patch_st: MagicMock, aspirin_data: dict
    ) -> None:
        """Profile page displays PubChem image when CID is available."""
        patch_st.session_state["selected_compound"] = aspirin_data
        from chemfuse.web.pages.profile import render

        render()
        patch_st.image.assert_called()
        call_args = str(patch_st.image.call_args)
        assert "2244" in call_args

    def test_render_shows_expandable_sections(
        self, patch_st: MagicMock, aspirin_data: dict
    ) -> None:
        """Profile page creates expander sections for additional data."""
        patch_st.session_state["selected_compound"] = aspirin_data
        from chemfuse.web.pages.profile import render

        render()
        assert patch_st.expander.called


class TestProfileIdentifiers:
    """Test identifier rendering."""

    def test_render_identifiers_shows_table(
        self, patch_st: MagicMock, aspirin_data: dict
    ) -> None:
        """Identifier section renders a dataframe."""
        from chemfuse.web.pages.profile import _render_identifiers

        _render_identifiers(aspirin_data)
        patch_st.dataframe.assert_called()

    def test_render_identifiers_includes_cid(
        self, patch_st: MagicMock, aspirin_data: dict
    ) -> None:
        """Identifier table includes the CID value."""
        from chemfuse.web.pages.profile import _render_identifiers

        _render_identifiers(aspirin_data)
        # CID 2244 should appear in data passed to dataframe
        call_args = patch_st.dataframe.call_args
        import pandas as pd
        df = call_args[0][0]
        assert isinstance(df, pd.DataFrame)
        assert "2244" in df["Value"].astype(str).values


class TestProfileProperties:
    """Test property rendering."""

    def test_render_properties_shows_table(
        self, patch_st: MagicMock, aspirin_data: dict
    ) -> None:
        """Properties section renders a table and chart."""
        from chemfuse.web.pages.profile import _render_properties

        with patch("chemfuse.web.components.property_chart.render_radar_chart"):
            _render_properties(aspirin_data)
        patch_st.dataframe.assert_called()


class TestProfileDruglikeness:
    """Test drug-likeness rendering."""

    def test_no_data_shows_caption(self, patch_st: MagicMock) -> None:
        """Drug-likeness section shows caption when no data available."""
        empty_compound: dict = {"smiles": ""}
        from chemfuse.web.pages.profile import _render_druglikeness

        _render_druglikeness(empty_compound)
        patch_st.caption.assert_called()

    def test_with_druglikeness_shows_badges(self, patch_st: MagicMock) -> None:
        """Drug-likeness section renders pass/fail badges."""
        compound = {
            "smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "druglikeness": {
                "filters": {
                    "Lipinski": {"pass": True},
                    "Veber": {"pass": True},
                }
            },
        }
        from chemfuse.web.pages.profile import _render_druglikeness

        _render_druglikeness(compound)
        patch_st.markdown.assert_called()


class TestProfileEnrichmentButtons:
    """Test enrich button rendering."""

    def test_enrich_buttons_displayed(
        self, patch_st: MagicMock, aspirin_data: dict
    ) -> None:
        """Profile page shows enrich buttons for each data source."""
        patch_st.session_state["selected_compound"] = aspirin_data
        from chemfuse.web.pages.profile import render

        render()
        # Buttons should be created for ChEMBL, BindingDB, SureChEMBL, OpenTargets
        all_calls = " ".join(str(c) for c in patch_st.button.call_args_list)
        assert "Enrich" in all_calls
