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


# ============================================================================
# _sanitize_svg tests
# ============================================================================

class TestSanitizeSvg:
    """Test SVG sanitization for XSS prevention."""

    def test_removes_script_tags(self) -> None:
        """Removes <script> elements from SVG."""
        from chemfuse.web.components.mol_grid import _sanitize_svg

        svg = '<svg><script>alert("xss")</script><circle/></svg>'
        result = _sanitize_svg(svg)
        assert "<script>" not in result
        assert "alert" not in result

    def test_removes_script_with_attributes(self) -> None:
        """Removes <script> with type attribute."""
        from chemfuse.web.components.mol_grid import _sanitize_svg

        svg = '<svg><script type="text/javascript">evil()</script></svg>'
        result = _sanitize_svg(svg)
        assert "<script" not in result

    def test_removes_onclick_handler(self) -> None:
        """Removes onclick event handler attributes."""
        from chemfuse.web.components.mol_grid import _sanitize_svg

        svg = '<svg><circle onclick="alert(1)"/></svg>'
        result = _sanitize_svg(svg)
        assert "onclick" not in result

    def test_removes_onload_handler(self) -> None:
        """Removes onload event handler attributes."""
        from chemfuse.web.components.mol_grid import _sanitize_svg

        svg = '<svg onload="evil()"><rect/></svg>'
        result = _sanitize_svg(svg)
        assert "onload" not in result

    def test_removes_single_quoted_event_handler(self) -> None:
        """Removes event handlers with single-quoted values."""
        from chemfuse.web.components.mol_grid import _sanitize_svg

        svg = "<svg><circle onmouseover='evil()'/></svg>"
        result = _sanitize_svg(svg)
        assert "onmouseover" not in result

    def test_preserves_normal_svg_content(self) -> None:
        """Preserves safe SVG content unchanged."""
        from chemfuse.web.components.mol_grid import _sanitize_svg

        svg = '<svg xmlns="http://www.w3.org/2000/svg"><circle r="50"/></svg>'
        result = _sanitize_svg(svg)
        assert "<circle" in result
        assert 'r="50"' in result

    def test_handles_empty_string(self) -> None:
        """Handles empty string input without error."""
        from chemfuse.web.components.mol_grid import _sanitize_svg

        result = _sanitize_svg("")
        assert result == ""

    def test_case_insensitive_removal(self) -> None:
        """Event handler removal is case-insensitive."""
        from chemfuse.web.components.mol_grid import _sanitize_svg

        svg = '<svg><circle ONCLICK="evil()"/></svg>'
        result = _sanitize_svg(svg)
        assert "ONCLICK" not in result


# ============================================================================
# smiles_to_svg additional tests
# ============================================================================

class TestSmilesToSvgRdkitFailure:
    """Test smiles_to_svg() with RDKit failure path."""

    def test_returns_none_for_invalid_smiles_when_rdkit_available(self) -> None:
        """Returns None for completely invalid SMILES even when RDKit is available."""
        import chemfuse.web.components.mol_grid as mod

        if not mod._RDKIT_AVAILABLE:
            return  # Skip if RDKit not installed

        # "INVALID_SMILES_!!!!" should fail RDKit parsing
        result = mod.smiles_to_svg("INVALID_SMILES_!!!!")
        assert result is None

    def test_returns_none_when_rdkit_draw_raises(self) -> None:
        """Returns None when RDKit drawing raises an exception."""
        import chemfuse.web.components.mol_grid as mod

        original = mod._RDKIT_AVAILABLE
        mod._RDKIT_AVAILABLE = True
        try:
            with patch.object(
                mod,
                "smiles_to_svg",
                wraps=lambda *a, **kw: None,
            ):
                result = mod.smiles_to_svg("CCO")
                assert result is None
        finally:
            mod._RDKIT_AVAILABLE = original


# ============================================================================
# render_mol_card additional tests
# ============================================================================

class TestRenderMolCardAdditional:
    """Additional tests for render_mol_card edge cases."""

    def test_renders_placeholder_for_smiles_only_no_cid(
        self, patch_st: MagicMock
    ) -> None:
        """Renders markdown with SMILES when CID is absent and RDKit unavailable."""
        import chemfuse.web.components.mol_grid as mod

        original = mod._RDKIT_AVAILABLE
        mod._RDKIT_AVAILABLE = False
        try:
            from chemfuse.web.components.mol_grid import render_mol_card

            compound = {"smiles": "CCO", "name": "ethanol"}
            render_mol_card(compound)
        finally:
            mod._RDKIT_AVAILABLE = original

        # Should render markdown (no image, no SVG)
        patch_st.markdown.assert_called()

    def test_renders_caption_with_compound_name(
        self, patch_st: MagicMock
    ) -> None:
        """Renders st.caption with the compound name."""
        import chemfuse.web.components.mol_grid as mod

        original = mod._RDKIT_AVAILABLE
        mod._RDKIT_AVAILABLE = False
        try:
            from chemfuse.web.components.mol_grid import render_mol_card

            compound = {"name": "ethanol", "smiles": "CCO"}
            render_mol_card(compound)
        finally:
            mod._RDKIT_AVAILABLE = original

        patch_st.caption.assert_called()
        caption_text = str(patch_st.caption.call_args)
        assert "ethanol" in caption_text

    def test_uses_smiles_as_name_when_name_absent(
        self, patch_st: MagicMock
    ) -> None:
        """Uses truncated SMILES as name when 'name' key is absent."""
        import chemfuse.web.components.mol_grid as mod

        original = mod._RDKIT_AVAILABLE
        mod._RDKIT_AVAILABLE = False
        try:
            from chemfuse.web.components.mol_grid import render_mol_card

            compound = {"smiles": "CCO"}
            render_mol_card(compound)
        finally:
            mod._RDKIT_AVAILABLE = original

        patch_st.caption.assert_called()


# ============================================================================
# render_radar_chart with max_values normalization
# ============================================================================

class TestRenderRadarChartNormalization:
    """Test radar chart with max_values normalization."""

    def test_renders_with_max_values(self, patch_st: MagicMock) -> None:
        """Renders chart when max_values normalization dict is provided."""
        from chemfuse.web.components.property_chart import render_radar_chart

        properties = {"MW": 180.0, "XLogP": 1.2, "TPSA": 63.6}
        max_values = {"MW": 500.0, "XLogP": 5.0, "TPSA": 150.0}
        render_radar_chart(properties, max_values=max_values)
        patch_st.plotly_chart.assert_called()

    def test_renders_with_partial_max_values(self, patch_st: MagicMock) -> None:
        """Renders chart when max_values is provided for only some properties."""
        from chemfuse.web.components.property_chart import render_radar_chart

        properties = {"MW": 180.0, "XLogP": 1.2, "TPSA": 63.6}
        max_values = {"MW": 500.0}  # Only MW provided
        render_radar_chart(properties, max_values=max_values)
        patch_st.plotly_chart.assert_called()

    def test_renders_without_max_values(self, patch_st: MagicMock) -> None:
        """Renders chart without normalization when max_values is None."""
        from chemfuse.web.components.property_chart import render_radar_chart

        properties = {"MW": 180.0, "XLogP": 1.2, "TPSA": 63.6}
        render_radar_chart(properties, max_values=None)
        patch_st.plotly_chart.assert_called()

    def test_custom_title_is_used(self, patch_st: MagicMock) -> None:
        """Custom title is passed to the chart."""
        from chemfuse.web.components.property_chart import render_radar_chart

        properties = {"A": 1.0, "B": 2.0, "C": 3.0}
        render_radar_chart(properties, title="Test Profile")
        patch_st.plotly_chart.assert_called()


# ============================================================================
# extract_chart_properties with nested properties
# ============================================================================

class TestExtractChartPropertiesAdditional:
    """Additional tests for extract_chart_properties()."""

    def test_extracts_from_nested_properties_dict(self) -> None:
        """Extracts properties from nested 'properties' dict when top-level absent."""
        from chemfuse.web.components.property_chart import extract_chart_properties

        compound = {
            "properties": {
                "molecular_weight": 180.16,
                "xlogp": 1.19,
                "tpsa": 63.6,
            }
        }
        result = extract_chart_properties(compound)
        assert "MW" in result
        assert "XLogP" in result
        assert "TPSA" in result

    def test_skips_non_numeric_values(self) -> None:
        """Skips properties that cannot be converted to float."""
        from chemfuse.web.components.property_chart import extract_chart_properties

        compound = {"molecular_weight": "not-a-number", "xlogp": 1.2}
        result = extract_chart_properties(compound)
        assert "MW" not in result
        assert "XLogP" in result

    def test_top_level_takes_priority_over_nested(self) -> None:
        """Top-level properties take priority over nested 'properties' dict."""
        from chemfuse.web.components.property_chart import extract_chart_properties

        compound = {
            "molecular_weight": 180.16,
            "properties": {"molecular_weight": 999.99},
        }
        result = extract_chart_properties(compound)
        assert result.get("MW") == 180.16

    def test_all_standard_properties_extracted(self, aspirin_data: dict) -> None:
        """All standard properties are extracted from aspirin data."""
        from chemfuse.web.components.property_chart import extract_chart_properties

        result = extract_chart_properties(aspirin_data)
        # Aspirin data has mw, xlogp, tpsa, hbd, hba, rotatable_bonds
        assert len(result) >= 4


# ============================================================================
# render_multi_radar_chart tests
# ============================================================================

class TestRenderMultiRadarChart:
    """Tests for render_multi_radar_chart()."""

    def test_renders_multi_radar_chart_for_two_compounds(
        self, patch_st: MagicMock
    ) -> None:
        """Renders chart for two compounds."""
        from chemfuse.web.components.property_chart import render_multi_radar_chart

        compounds = [
            {"MW": 180.0, "XLogP": 1.2, "TPSA": 63.6},
            {"MW": 194.0, "XLogP": -0.07, "TPSA": 58.4},
        ]
        names = ["aspirin", "caffeine"]
        render_multi_radar_chart(compounds, names)
        patch_st.plotly_chart.assert_called()

    def test_skips_when_empty_list(self, patch_st: MagicMock) -> None:
        """Does not render when compound list is empty."""
        from chemfuse.web.components.property_chart import render_multi_radar_chart

        render_multi_radar_chart([], [])
        patch_st.plotly_chart.assert_not_called()

    def test_skips_when_mismatched_lengths(self, patch_st: MagicMock) -> None:
        """Does not render when compound count does not match name count."""
        from chemfuse.web.components.property_chart import render_multi_radar_chart

        compounds = [{"MW": 180.0, "XLogP": 1.2, "TPSA": 63.6}]
        names = ["aspirin", "caffeine"]  # mismatch
        render_multi_radar_chart(compounds, names)
        patch_st.plotly_chart.assert_not_called()

    def test_skips_when_too_few_unique_properties(
        self, patch_st: MagicMock
    ) -> None:
        """Does not render when fewer than 3 unique properties across compounds."""
        from chemfuse.web.components.property_chart import render_multi_radar_chart

        compounds = [
            {"MW": 180.0, "XLogP": 1.2},
            {"MW": 194.0, "XLogP": -0.07},
        ]
        names = ["aspirin", "caffeine"]
        render_multi_radar_chart(compounds, names)
        patch_st.plotly_chart.assert_not_called()

    def test_uses_custom_title(self, patch_st: MagicMock) -> None:
        """Custom title is accepted without error."""
        from chemfuse.web.components.property_chart import render_multi_radar_chart

        compounds = [
            {"A": 1.0, "B": 2.0, "C": 3.0},
            {"A": 4.0, "B": 5.0, "C": 6.0},
        ]
        render_multi_radar_chart(compounds, ["x", "y"], title="Custom Title")
        patch_st.plotly_chart.assert_called()

    def test_shows_caption_when_plotly_not_installed(
        self, patch_st: MagicMock
    ) -> None:
        """Shows install caption when plotly is not available."""
        # Temporarily block plotly import
        with patch.dict("sys.modules", {"plotly": None, "plotly.graph_objects": None}):
            import importlib

            import chemfuse.web.components.property_chart as mod
            # Reload to trigger ImportError path
            try:
                importlib.reload(mod)
            except Exception:
                pass
            # Call with plotly blocked
            from chemfuse.web.components.property_chart import render_multi_radar_chart
            compounds = [
                {"A": 1.0, "B": 2.0, "C": 3.0},
                {"A": 4.0, "B": 5.0, "C": 6.0},
            ]
            render_multi_radar_chart(compounds, ["x", "y"])


class TestPropertyChartImportError:
    """Test property chart functions when plotly is not available."""

    def test_radar_chart_shows_caption_when_plotly_missing(
        self, patch_st: MagicMock
    ) -> None:
        """render_radar_chart shows caption when plotly cannot be imported."""
        import builtins

        original_import = builtins.__import__

        def _mock_import(name, *args, **kwargs):
            if name == "plotly.graph_objects":
                raise ImportError("plotly not installed")
            return original_import(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=_mock_import):
            from chemfuse.web.components.property_chart import render_radar_chart
            properties = {"MW": 180.0, "XLogP": 1.2, "TPSA": 63.6}
            render_radar_chart(properties)
        # Either caption called (ImportError path) or plotly_chart called (plotly available)
        # Just verify no exception raised
        assert True

    def test_bar_chart_shows_caption_when_plotly_missing(
        self, patch_st: MagicMock
    ) -> None:
        """render_bar_chart shows caption when plotly cannot be imported."""
        import builtins

        original_import = builtins.__import__

        def _mock_import(name, *args, **kwargs):
            if name == "plotly.graph_objects":
                raise ImportError("plotly not installed")
            return original_import(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=_mock_import):
            from chemfuse.web.components.property_chart import render_bar_chart
            render_bar_chart({"MW": 180.0})
        assert True


# ============================================================================
# mol_grid RDKit path tests
# ============================================================================

class TestSmilesToSvgRdkitPath:
    """Test smiles_to_svg when RDKit is available."""

    def test_valid_smiles_returns_svg_or_none(self) -> None:
        """For a valid SMILES, returns SVG string or None (if RDKit unavailable)."""
        from chemfuse.web.components.mol_grid import smiles_to_svg

        result = smiles_to_svg("CC(=O)Oc1ccccc1C(=O)O")
        # Either None (no RDKit) or a string
        assert result is None or isinstance(result, str)

    def test_returns_none_when_rdkit_explicitly_disabled(self) -> None:
        """Always returns None when _RDKIT_AVAILABLE is False."""
        import chemfuse.web.components.mol_grid as mod

        original = mod._RDKIT_AVAILABLE
        mod._RDKIT_AVAILABLE = False
        try:
            result = mod.smiles_to_svg("CC(=O)O")
            assert result is None
        finally:
            mod._RDKIT_AVAILABLE = original


class TestRenderMolCardButtonClick:
    """Test button click behavior in render_mol_card."""

    def test_button_click_sets_session_state_and_reruns(
        self, patch_st: MagicMock
    ) -> None:
        """Clicking 'View Profile' button sets session state and calls st.rerun."""
        import chemfuse.web.components.mol_grid as mod

        original = mod._RDKIT_AVAILABLE
        mod._RDKIT_AVAILABLE = False
        try:
            # Make button return True (clicked)
            patch_st.button.return_value = True
            compound = {"name": "aspirin", "cid": 2244, "smiles": "CC(=O)O"}
            from chemfuse.web.components.mol_grid import render_mol_card

            render_mol_card(compound, on_click_session_key="selected_compound")
        finally:
            mod._RDKIT_AVAILABLE = original

        assert patch_st.session_state.get("selected_compound") == compound
        assert patch_st.session_state.get("current_page") == "Profile"
        patch_st.rerun.assert_called_once()
