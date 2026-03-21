"""Tests for the Compare page."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest


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

    def test_render_with_too_few_selected_shows_info(
        self, patch_st: MagicMock, compound_list: list[dict]
    ) -> None:
        """Shows info message when fewer than 2 compounds selected."""
        patch_st.session_state["session_compounds"] = compound_list
        patch_st.multiselect.return_value = [0]  # Only 1 selected
        from chemfuse.web.pages.compare import render

        render()
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
