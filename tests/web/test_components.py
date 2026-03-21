"""Tests for web UI components: mol_grid, property_chart, filters."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

# ============================================================================
# mol_grid component tests
# ============================================================================

class TestMolGridIsRdkitAvailable:
    """Test RDKit availability check."""

    def test_returns_bool(self) -> None:
        from chemfuse.web.components.mol_grid import is_rdkit_available

        result = is_rdkit_available()
        assert isinstance(result, bool)


class TestSmilesToSvg:
    """Test SVG generation from SMILES."""

    def test_returns_none_when_no_rdkit(self) -> None:
        """Returns None when RDKit is not available."""
        with patch("chemfuse.web.components.mol_grid._RDKIT_AVAILABLE", False):
            from chemfuse.web.components import mol_grid
            old_val = mol_grid._RDKIT_AVAILABLE
            mol_grid._RDKIT_AVAILABLE = False
            result = mol_grid.smiles_to_svg("CCO")
            mol_grid._RDKIT_AVAILABLE = old_val
        assert result is None

    def test_returns_none_for_empty_smiles(self) -> None:
        """Returns None for empty SMILES string."""
        from chemfuse.web.components.mol_grid import smiles_to_svg

        result = smiles_to_svg("")
        assert result is None

    def test_returns_none_when_rdkit_unavailable_at_module_level(self) -> None:
        """Returns None when RDKit is marked unavailable."""
        import chemfuse.web.components.mol_grid as mod

        original = mod._RDKIT_AVAILABLE
        mod._RDKIT_AVAILABLE = False
        try:
            result = mod.smiles_to_svg("CCO")
            assert result is None
        finally:
            mod._RDKIT_AVAILABLE = original


class TestGetPubchemPngUrl:
    """Test PubChem PNG URL generation."""

    def test_returns_url_for_valid_cid(self) -> None:
        """Returns valid PubChem PNG URL for a CID."""
        from chemfuse.web.components.mol_grid import get_pubchem_png_url

        url = get_pubchem_png_url(2244)
        assert url is not None
        assert "2244" in url
        assert "pubchem" in url.lower()
        assert url.endswith("PNG")

    def test_returns_none_for_none_cid(self) -> None:
        """Returns None when CID is None."""
        from chemfuse.web.components.mol_grid import get_pubchem_png_url

        assert get_pubchem_png_url(None) is None


class TestRenderMolCard:
    """Test single molecule card rendering."""

    def test_renders_pubchem_image_when_cid_available(
        self, patch_st: MagicMock, aspirin_data: dict
    ) -> None:
        """Renders PubChem PNG image when CID is available and RDKit is absent."""
        import chemfuse.web.components.mol_grid as mod

        original = mod._RDKIT_AVAILABLE
        mod._RDKIT_AVAILABLE = False
        try:
            mod.render_mol_card(aspirin_data)
        finally:
            mod._RDKIT_AVAILABLE = original

        patch_st.image.assert_called()
        url_arg = str(patch_st.image.call_args)
        assert "2244" in url_arg

    def test_renders_placeholder_when_no_smiles_and_no_cid(
        self, patch_st: MagicMock
    ) -> None:
        """Renders placeholder when no structure info available."""
        from chemfuse.web.components.mol_grid import render_mol_card

        compound: dict = {"name": "unknown_compound"}
        render_mol_card(compound)
        patch_st.markdown.assert_called()

    def test_renders_select_button_when_key_provided(
        self, patch_st: MagicMock, aspirin_data: dict
    ) -> None:
        """Renders 'View Profile' button when on_click_session_key is provided."""
        import chemfuse.web.components.mol_grid as mod

        mod._RDKIT_AVAILABLE = False
        from chemfuse.web.components.mol_grid import render_mol_card

        render_mol_card(aspirin_data, on_click_session_key="selected_compound")
        button_calls = " ".join(str(c) for c in patch_st.button.call_args_list)
        assert "View Profile" in button_calls


class TestRenderMolGrid:
    """Test molecule grid rendering."""

    def test_empty_list_shows_info(self, patch_st: MagicMock) -> None:
        """Empty compound list shows info message."""
        from chemfuse.web.components.mol_grid import render_mol_grid

        render_mol_grid([])
        patch_st.info.assert_called()

    def test_renders_columns_for_compounds(
        self, patch_st: MagicMock, compound_list: list[dict]
    ) -> None:
        """Creates st.columns for molecule grid layout."""
        from chemfuse.web.components.mol_grid import render_mol_grid

        render_mol_grid(compound_list)
        patch_st.columns.assert_called()

    def test_max_compounds_limit(self, patch_st: MagicMock) -> None:
        """Grid respects max_compounds limit and shows caption."""
        many_compounds = [{"smiles": f"C{i}", "name": f"c{i}"} for i in range(100)]
        from chemfuse.web.components.mol_grid import render_mol_grid

        render_mol_grid(many_compounds, max_compounds=10)
        patch_st.caption.assert_called()
        # Caption should mention truncation
        caption_text = " ".join(str(c) for c in patch_st.caption.call_args_list)
        assert "10" in caption_text


# ============================================================================
# property_chart component tests
# ============================================================================

class TestRenderRadarChart:
    """Test radar chart rendering."""

    def test_renders_plotly_chart(self, patch_st: MagicMock) -> None:
        """Renders a Plotly radar chart for valid property data."""
        from chemfuse.web.components.property_chart import render_radar_chart

        properties = {"MW": 180.0, "XLogP": 1.2, "TPSA": 63.6, "HBD": 1.0, "HBA": 4.0}
        render_radar_chart(properties)
        patch_st.plotly_chart.assert_called()

    def test_skips_chart_with_too_few_properties(self, patch_st: MagicMock) -> None:
        """Shows caption and skips chart when fewer than 3 properties."""
        from chemfuse.web.components.property_chart import render_radar_chart

        render_radar_chart({"MW": 180.0, "XLogP": 1.2})
        patch_st.plotly_chart.assert_not_called()
        patch_st.caption.assert_called()


class TestRenderBarChart:
    """Test bar chart rendering."""

    def test_renders_bar_chart(self, patch_st: MagicMock) -> None:
        """Renders a Plotly bar chart for property data."""
        from chemfuse.web.components.property_chart import render_bar_chart

        render_bar_chart({"MW": 180.0, "XLogP": 1.2})
        patch_st.plotly_chart.assert_called()

    def test_skips_empty_data(self, patch_st: MagicMock) -> None:
        """Skips rendering for empty property dict."""
        from chemfuse.web.components.property_chart import render_bar_chart

        render_bar_chart({})
        patch_st.plotly_chart.assert_not_called()


class TestExtractChartProperties:
    """Test property extraction for charts."""

    def test_extracts_numeric_properties(self, aspirin_data: dict) -> None:
        """Extracts numeric properties from compound dict."""
        from chemfuse.web.components.property_chart import extract_chart_properties

        result = extract_chart_properties(aspirin_data)
        assert "MW" in result
        assert isinstance(result["MW"], float)

    def test_returns_empty_for_empty_compound(self) -> None:
        """Returns empty dict for empty compound."""
        from chemfuse.web.components.property_chart import extract_chart_properties

        result = extract_chart_properties({})
        assert result == {}

    def test_extracts_from_nested_properties(self) -> None:
        """Extracts properties from nested 'properties' key."""
        from chemfuse.web.components.property_chart import extract_chart_properties

        compound = {"properties": {"molecular_weight": 180.16, "xlogp": 1.19}}
        result = extract_chart_properties(compound)
        assert "MW" in result or "XLogP" in result


# ============================================================================
# filters component tests
# ============================================================================

class TestRenderPropertyFilters:
    """Test property range filter widgets."""

    def test_returns_dict_of_ranges(self, patch_st: MagicMock) -> None:
        """Returns a dict of property ranges."""
        from chemfuse.web.components.filters import render_property_filters

        result = render_property_filters()
        assert isinstance(result, dict)

    def test_mw_filter_included_by_default(self, patch_st: MagicMock) -> None:
        """MW filter is included when show_mw=True."""
        patch_st.slider.return_value = (100.0, 500.0)
        from chemfuse.web.components.filters import render_property_filters

        result = render_property_filters(show_mw=True)
        assert "molecular_weight" in result

    def test_mw_filter_excluded(self, patch_st: MagicMock) -> None:
        """MW filter is excluded when show_mw=False."""
        from chemfuse.web.components.filters import render_property_filters

        result = render_property_filters(show_mw=False, show_logp=False, show_tpsa=False)
        assert "molecular_weight" not in result

    def test_range_values_come_from_slider(self, patch_st: MagicMock) -> None:
        """Filter values come from st.slider return value."""
        patch_st.slider.return_value = (50.0, 300.0)
        from chemfuse.web.components.filters import render_property_filters

        result = render_property_filters(show_mw=True, show_logp=False, show_tpsa=False)
        assert result["molecular_weight"] == (50.0, 300.0)


class TestRenderDruglikenessCheckboxes:
    """Test drug-likeness filter checkboxes."""

    def test_returns_list_of_selected_filters(self, patch_st: MagicMock) -> None:
        """Returns a list of selected filter names."""
        patch_st.checkbox.return_value = True
        from chemfuse.web.components.filters import render_druglikeness_checkboxes

        result = render_druglikeness_checkboxes()
        assert isinstance(result, list)
        assert len(result) > 0

    def test_returns_empty_list_when_none_selected(self, patch_st: MagicMock) -> None:
        """Returns empty list when no checkboxes are checked."""
        patch_st.checkbox.return_value = False
        from chemfuse.web.components.filters import render_druglikeness_checkboxes

        result = render_druglikeness_checkboxes()
        assert result == []
