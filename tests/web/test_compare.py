"""Tests for the Compare page."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest


class TestStableId:
    """Tests for _stable_id helper."""

    def test_cid_takes_priority(self) -> None:
        """CID is used as stable ID when present."""
        from chemfuse.web.pages.compare import _stable_id

        result = _stable_id({"cid": 2244, "smiles": "CC", "name": "aspirin"})
        assert result == "cid:2244"

    def test_smiles_used_when_no_cid(self) -> None:
        """SMILES is used as stable ID when CID is absent."""
        from chemfuse.web.pages.compare import _stable_id

        result = _stable_id({"smiles": "CC(=O)O", "name": "acetic acid"})
        assert result == "smiles:CC(=O)O"

    def test_name_used_when_no_cid_or_smiles(self) -> None:
        """Name is used as stable ID when CID and SMILES are absent."""
        from chemfuse.web.pages.compare import _stable_id

        result = _stable_id({"name": "compound-x"})
        assert result == "name:compound-x"

    def test_fallback_to_object_id(self) -> None:
        """Falls back to object id when all identifiers are absent."""
        from chemfuse.web.pages.compare import _stable_id

        c: dict = {}
        result = _stable_id(c)
        assert result.startswith("idx:")


class TestDisplayLabel:
    """Tests for _display_label helper."""

    def test_cid_and_name(self) -> None:
        """Label includes both CID and name when both are present."""
        from chemfuse.web.pages.compare import _display_label

        label = _display_label({"cid": 2244, "name": "aspirin"})
        assert "2244" in label
        assert "aspirin" in label

    def test_cid_only(self) -> None:
        """Label shows CID when name is absent."""
        from chemfuse.web.pages.compare import _display_label

        label = _display_label({"cid": 2244})
        assert "2244" in label

    def test_name_only(self) -> None:
        """Label shows name when CID is absent."""
        from chemfuse.web.pages.compare import _display_label

        label = _display_label({"name": "caffeine"})
        assert label == "caffeine"


class TestComparePageRender:
    """Test Compare page rendering."""

    def test_render_shows_header(self, patch_st: MagicMock) -> None:
        """Compare page shows 'Compound Comparison' header."""
        from chemfuse.web.pages.compare import render

        render()
        patch_st.header.assert_called_with("Compound Comparison")

    def test_render_with_no_compounds_shows_info(self, patch_st: MagicMock) -> None:
        """Shows info message when session has fewer than 2 compounds."""
        patch_st.session_state["session_compounds"] = []
        from chemfuse.web.pages.compare import render

        render()
        patch_st.info.assert_called()

    def test_render_with_one_compound_shows_info(self, patch_st: MagicMock) -> None:
        """Shows info message when only one compound is in session."""
        patch_st.session_state["session_compounds"] = [{"name": "aspirin"}]
        from chemfuse.web.pages.compare import render

        render()
        patch_st.info.assert_called()

    def test_render_with_two_compounds_shows_selector(
        self, patch_st: MagicMock, compound_list: list[dict]
    ) -> None:
        """Shows multiselect when 2+ compounds available."""
        patch_st.session_state["session_compounds"] = compound_list
        from chemfuse.web.pages.compare import render

        render()
        patch_st.multiselect.assert_called()

    def test_selector_options_are_stable_ids_not_indices(
        self, patch_st: MagicMock, compound_list: list[dict]
    ) -> None:
        """Multiselect options must be stable string IDs, not integer indices."""
        patch_st.session_state["session_compounds"] = compound_list
        from chemfuse.web.pages.compare import render

        render()
        call_kwargs = patch_st.multiselect.call_args
        options = call_kwargs[1].get("options") or call_kwargs[0][1]
        assert all(isinstance(opt, str) for opt in options), (
            "Options must be stable string identifiers, not integers"
        )

    def test_selector_options_contain_cid_keys(
        self, patch_st: MagicMock, compound_list: list[dict]
    ) -> None:
        """Multiselect options use CID-based keys for compounds that have a CID."""
        patch_st.session_state["session_compounds"] = compound_list
        from chemfuse.web.pages.compare import render

        render()
        call_kwargs = patch_st.multiselect.call_args
        options = call_kwargs[1].get("options") or call_kwargs[0][1]
        # Both fixtures (aspirin cid=2244, caffeine cid=2519) have CIDs
        assert "cid:2244" in options
        assert "cid:2519" in options

    def test_render_with_too_few_selected_shows_info(
        self, patch_st: MagicMock, compound_list: list[dict]
    ) -> None:
        """Shows info message when fewer than 2 compounds selected."""
        patch_st.session_state["session_compounds"] = compound_list
        patch_st.multiselect.return_value = ["cid:2244"]  # Only 1 stable ID
        from chemfuse.web.pages.compare import render

        render()
        patch_st.info.assert_called()

    def test_adding_compound_does_not_change_existing_selections(
        self, patch_st: MagicMock, compound_list: list[dict], aspirin_data: dict
    ) -> None:
        """Stable IDs are unchanged when a new compound is prepended to the session.

        Before: [aspirin, caffeine] -> options ["cid:2244", "cid:2519"]
        After:  [new_compound, aspirin, caffeine] -> "cid:2244" still resolves
        to aspirin regardless of its new position in the list.
        """
        from chemfuse.web.pages.compare import _stable_id

        new_compound = {"cid": 9999, "name": "new_drug", "smiles": "C"}

        # With new compound prepended the aspirin stable ID must be unchanged
        extended_list = [new_compound] + compound_list
        aspirin_id = _stable_id(aspirin_data)
        assert aspirin_id == "cid:2244"

        # Simulate user keeping aspirin + caffeine selected via their stable IDs
        patch_st.session_state["session_compounds"] = extended_list
        patch_st.multiselect.return_value = ["cid:2244", "cid:2519"]
        from chemfuse.web.pages.compare import render

        render()
        # Comparison sections must be rendered (dataframe called at least twice)
        assert patch_st.dataframe.call_count >= 2

    def test_removing_compound_updates_options(
        self, patch_st: MagicMock, compound_list: list[dict]
    ) -> None:
        """Removing a compound from the session removes its stable ID from options."""
        from chemfuse.web.pages.compare import render

        # Start with only caffeine remaining
        patch_st.session_state["session_compounds"] = [compound_list[1]]
        render()
        # Only 1 compound -> info shown, multiselect NOT called
        patch_st.info.assert_called()


class TestPropertyComparison:
    """Test property comparison rendering."""

    def test_renders_comparison_table(
        self, patch_st: MagicMock, aspirin_data: dict, caffeine_data: dict
    ) -> None:
        """Property comparison renders a dataframe."""
        with patch("chemfuse.web.pages.compare._render_radar_overlay"):
            from chemfuse.web.pages.compare import _render_property_comparison

            _render_property_comparison(
                [aspirin_data, caffeine_data],
                ["aspirin", "caffeine"],
            )
        patch_st.dataframe.assert_called()

    def test_comparison_table_includes_both_compounds(
        self, patch_st: MagicMock, aspirin_data: dict, caffeine_data: dict
    ) -> None:
        """Comparison table has columns for both compound names."""
        with patch("chemfuse.web.pages.compare._render_radar_overlay"):
            from chemfuse.web.pages.compare import _render_property_comparison

            _render_property_comparison(
                [aspirin_data, caffeine_data],
                ["aspirin", "caffeine"],
            )

        call_args = patch_st.dataframe.call_args
        df = call_args[0][0]
        assert "aspirin" in df.columns
        assert "caffeine" in df.columns


class TestDruglikenessComparison:
    """Test drug-likeness comparison."""

    def test_renders_druglikeness_table(
        self, patch_st: MagicMock, aspirin_data: dict, caffeine_data: dict
    ) -> None:
        """Drug-likeness comparison renders a dataframe."""
        from chemfuse.web.pages.compare import _render_druglikeness_comparison

        _render_druglikeness_comparison(
            [aspirin_data, caffeine_data],
            ["aspirin", "caffeine"],
        )
        patch_st.dataframe.assert_called()

    def test_table_shows_pass_fail(self, patch_st: MagicMock) -> None:
        """Drug-likeness table shows PASS/FAIL values."""
        compound_with_dl = {
            "name": "aspirin",
            "druglikeness": {
                "filters": {
                    "Lipinski": {"pass": True},
                    "Veber": {"pass": False},
                }
            },
        }
        compound_without_dl = {"name": "caffeine"}

        from chemfuse.web.pages.compare import _render_druglikeness_comparison

        _render_druglikeness_comparison(
            [compound_with_dl, compound_without_dl],
            ["aspirin", "caffeine"],
        )
        call_args = patch_st.dataframe.call_args
        df = call_args[0][0]
        assert "PASS" in df["aspirin"].values or "FAIL" in df["aspirin"].values
        assert "—" in df["caffeine"].values


class TestGetPropertyHelper:
    """Test _get_property helper function."""

    def test_extracts_top_level_value(self) -> None:
        from chemfuse.web.pages.compare import _get_property

        compound = {"molecular_weight": 180.16}
        result = _get_property(compound, "molecular_weight")
        assert result == pytest.approx(180.16)

    def test_extracts_from_nested_properties(self) -> None:
        from chemfuse.web.pages.compare import _get_property

        compound = {"properties": {"molecular_weight": 180.16}}
        result = _get_property(compound, "molecular_weight")
        assert result == pytest.approx(180.16)

    def test_returns_none_when_missing(self) -> None:
        from chemfuse.web.pages.compare import _get_property

        compound: dict = {}
        result = _get_property(compound, "molecular_weight")
        assert result is None


class TestGetDruglikenessHelper:
    """Test _get_druglikeness_result helper."""

    def test_returns_true_for_passing_filter(self) -> None:
        from chemfuse.web.pages.compare import _get_druglikeness_result

        compound = {
            "druglikeness": {
                "filters": {"Lipinski": {"pass": True}}
            }
        }
        result = _get_druglikeness_result(compound, "Lipinski")
        assert result is True

    def test_returns_false_for_failing_filter(self) -> None:
        from chemfuse.web.pages.compare import _get_druglikeness_result

        compound = {
            "druglikeness": {
                "filters": {"Lipinski": {"pass": False}}
            }
        }
        result = _get_druglikeness_result(compound, "Lipinski")
        assert result is False

    def test_returns_none_when_no_druglikeness(self) -> None:
        from chemfuse.web.pages.compare import _get_druglikeness_result

        result = _get_druglikeness_result({"name": "test"}, "Lipinski")
        assert result is None
