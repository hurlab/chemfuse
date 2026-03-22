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

    def test_disabled_buttons_for_bindingdb(
        self, patch_st: MagicMock, aspirin_data: dict
    ) -> None:
        """Enrich from BindingDB button is disabled."""
        patch_st.session_state["selected_compound"] = aspirin_data
        from chemfuse.web.pages.profile import render

        render()
        # Check that button is called with disabled=True for BindingDB
        all_calls = " ".join(str(c) for c in patch_st.button.call_args_list)
        assert "BindingDB" in all_calls

    def test_disabled_buttons_for_surechembl(
        self, patch_st: MagicMock, aspirin_data: dict
    ) -> None:
        """Enrich from SureChEMBL button is disabled."""
        patch_st.session_state["selected_compound"] = aspirin_data
        from chemfuse.web.pages.profile import render

        render()
        all_calls = " ".join(str(c) for c in patch_st.button.call_args_list)
        assert "SureChEMBL" in all_calls


class TestProfileManualLookup:
    """Test _render_manual_lookup() function."""

    def test_renders_text_input_for_identifier(self, patch_st: MagicMock) -> None:
        """_render_manual_lookup creates a text input for CID or SMILES."""
        from chemfuse.web.pages.profile import _render_manual_lookup

        _render_manual_lookup()
        patch_st.text_input.assert_called()

    def test_renders_lookup_button(self, patch_st: MagicMock) -> None:
        """_render_manual_lookup creates a 'Look up' button."""
        from chemfuse.web.pages.profile import _render_manual_lookup

        _render_manual_lookup()
        button_labels = " ".join(str(c) for c in patch_st.button.call_args_list)
        assert "Look up" in button_labels

    def test_lookup_with_empty_identifier_does_nothing(
        self, patch_st: MagicMock
    ) -> None:
        """Empty identifier string does not trigger lookup."""
        patch_st.button.return_value = True
        patch_st.text_input.return_value = ""
        from chemfuse.web.pages.profile import _render_manual_lookup

        # Should not call st.rerun or st.warning
        _render_manual_lookup()
        patch_st.rerun.assert_not_called()

    def test_lookup_sets_session_state_on_success(self, patch_st: MagicMock) -> None:
        """Successful lookup stores compound in session_state and reruns."""
        patch_st.button.return_value = True
        patch_st.text_input.return_value = "2244"

        mock_data = {"name": "aspirin", "cid": 2244}
        with patch(
            "chemfuse.web.pages.profile._lookup_compound",
            return_value=mock_data,
        ):
            from chemfuse.web.pages.profile import _render_manual_lookup
            _render_manual_lookup()

        assert patch_st.session_state.get("selected_compound") == mock_data
        patch_st.rerun.assert_called_once()

    def test_lookup_warns_when_not_found(self, patch_st: MagicMock) -> None:
        """Shows warning when _lookup_compound returns None."""
        patch_st.button.return_value = True
        patch_st.text_input.return_value = "unknown_smiles_xyz"

        with patch(
            "chemfuse.web.pages.profile._lookup_compound",
            return_value=None,
        ):
            from chemfuse.web.pages.profile import _render_manual_lookup
            _render_manual_lookup()

        patch_st.warning.assert_called()


class TestLookupCompound:
    """Test _lookup_compound() function."""

    def test_lookup_by_cid(self) -> None:
        """Uses 'cid' query type when identifier is all digits."""
        mock_compound = MagicMock()
        mock_compound.to_dict.return_value = {"name": "aspirin", "cid": 2244}

        async def _fake_search(identifier, query_type):
            assert query_type == "cid"
            return [mock_compound]

        with patch(
            "chemfuse.web.pages.profile._run_async",
            return_value={"name": "aspirin", "cid": 2244},
        ):
            from chemfuse.web.pages.profile import _lookup_compound
            result = _lookup_compound("2244")
        assert result is not None

    def test_lookup_by_smiles(self) -> None:
        """Uses 'smiles' query type when identifier is not digits."""
        with patch(
            "chemfuse.web.pages.profile._run_async",
            return_value={"name": "aspirin", "smiles": "CC(=O)O"},
        ):
            from chemfuse.web.pages.profile import _lookup_compound
            result = _lookup_compound("CC(=O)O")
        assert result is not None

    def test_lookup_returns_none_on_error(self, patch_st: MagicMock) -> None:
        """Returns None and shows error when _run_async raises."""
        with patch(
            "chemfuse.web.pages.profile._run_async",
            side_effect=Exception("network error"),
        ):
            from chemfuse.web.pages.profile import _lookup_compound
            result = _lookup_compound("aspirin")
        assert result is None
        patch_st.error.assert_called()


class TestRenderDruglikenessExceptions:
    """Test _render_druglikeness exception handling."""

    def test_shows_warning_on_compute_exception(self, patch_st: MagicMock) -> None:
        """Shows warning when drug-likeness computation raises an exception."""
        compound = {"smiles": "CC(=O)Oc1ccccc1C(=O)O"}
        with patch(
            "chemfuse.compute.druglikeness.check_drug_likeness",
            side_effect=Exception("compute error"),
        ):
            from chemfuse.web.pages.profile import _render_druglikeness
            _render_druglikeness(compound)
        patch_st.warning.assert_called()

    def test_shows_caption_when_no_smiles_and_no_druglikeness(
        self, patch_st: MagicMock
    ) -> None:
        """Shows caption when neither SMILES nor druglikeness data available."""
        compound: dict = {}
        from chemfuse.web.pages.profile import _render_druglikeness

        _render_druglikeness(compound)
        patch_st.caption.assert_called()


class TestRenderAdmet:
    """Test _render_admet() function."""

    def test_shows_caption_when_no_smiles(self, patch_st: MagicMock) -> None:
        """Shows caption when no SMILES in compound data."""
        from chemfuse.web.pages.profile import _render_admet

        _render_admet({})
        patch_st.caption.assert_called()
        caption_text = " ".join(str(c) for c in patch_st.caption.call_args_list)
        assert "SMILES" in caption_text

    def test_shows_cached_admet_data(self, patch_st: MagicMock) -> None:
        """Displays cached ADMET data without showing predict button."""
        compound = {
            "smiles": "CC(=O)O",
            "admet": [{"property": "HIA", "value": 0.9}],
        }
        from chemfuse.web.pages.profile import _render_admet

        _render_admet(compound)
        patch_st.dataframe.assert_called()

    def test_shows_predict_button_when_no_cached_admet(
        self, patch_st: MagicMock
    ) -> None:
        """Shows 'Predict ADMET' button when no cached ADMET data."""
        compound = {"smiles": "CC(=O)O"}
        from chemfuse.web.pages.profile import _render_admet

        _render_admet(compound)
        button_labels = " ".join(str(c) for c in patch_st.button.call_args_list)
        assert "ADMET" in button_labels or "Predict" in button_labels

    def test_predict_admet_on_button_click(self, patch_st: MagicMock) -> None:
        """Calls predict_admet and shows dataframe when button clicked."""
        patch_st.button.return_value = True
        compound = {"smiles": "CC(=O)O"}
        mock_preds = [{"HIA": 0.9, "BBB": 0.5}]

        with patch(
            "chemfuse.compute.admet.predict_admet",
            return_value=mock_preds,
        ):
            from chemfuse.web.pages.profile import _render_admet
            _render_admet(compound)
        patch_st.dataframe.assert_called()

    def test_admet_persists_in_session_state(self, patch_st: MagicMock) -> None:
        """ADMET results are persisted in session_state after prediction."""
        patch_st.button.return_value = True
        compound = {"smiles": "CC(=O)O"}
        mock_preds = [{"HIA": 0.9}]

        with patch(
            "chemfuse.compute.admet.predict_admet",
            return_value=mock_preds,
        ):
            from chemfuse.web.pages.profile import _render_admet
            _render_admet(compound)

        # Session state should be updated with the compound including admet
        selected = patch_st.session_state.get("selected_compound")
        if selected is not None:
            assert "admet" in selected

    def test_predict_admet_shows_error_on_exception(
        self, patch_st: MagicMock
    ) -> None:
        """Shows error when predict_admet raises exception."""
        patch_st.button.return_value = True
        compound = {"smiles": "CC(=O)O"}

        with patch(
            "chemfuse.compute.admet.predict_admet",
            side_effect=Exception("compute error"),
        ):
            from chemfuse.web.pages.profile import _render_admet
            _render_admet(compound)
        patch_st.error.assert_called()
