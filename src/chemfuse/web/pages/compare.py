"""Compare page - side-by-side compound comparison."""

from __future__ import annotations

from typing import Any

import pandas as pd
import streamlit as st

_PROPERTIES_TO_COMPARE = [
    ("molecular_weight", "MW (Da)"),
    ("xlogp", "XLogP"),
    ("tpsa", "TPSA (A²)"),
    ("hbd_count", "HBD"),
    ("hba_count", "HBA"),
    ("rotatable_bonds", "Rotatable Bonds"),
    ("heavy_atom_count", "Heavy Atoms"),
    ("complexity", "Complexity"),
]

_DRUGLIKENESS_FILTERS = ["Lipinski", "Veber", "Ghose", "Egan", "Muegge"]


def render() -> None:
    """Render the Compare page."""
    st.header("Compound Comparison")
    st.caption("Select 2–4 compounds from your session for side-by-side comparison.")

    session_compounds = st.session_state.get("session_compounds", [])
    if len(session_compounds) < 2:
        st.info(
            "Not enough compounds in session. "
            "Run a search to populate the compound pool (need at least 2)."
        )
        return

    # Compound selection
    compound_names = [
        c.get("name") or c.get("smiles", "Unknown")[:20]
        for c in session_compounds
    ]
    selected_indices = st.multiselect(
        "Select compounds to compare (2–4)",
        options=list(range(len(session_compounds))),
        format_func=lambda i: compound_names[i],
        key="compare_selection",
        max_selections=4,
    )

    if len(selected_indices) < 2:
        st.info("Select at least 2 compounds.")
        return

    selected = [session_compounds[i] for i in selected_indices]
    names = [compound_names[i] for i in selected_indices]

    # Render comparison sections
    _render_property_comparison(selected, names)
    _render_druglikeness_comparison(selected, names)
    _render_admet_comparison(selected, names)


def _render_property_comparison(
    compounds: list[dict[str, Any]],
    names: list[str],
) -> None:
    """Render side-by-side property comparison table and radar chart."""
    st.subheader("Property Comparison")

    rows: list[dict[str, Any]] = []
    chart_data: dict[str, list[float | None]] = {}

    for prop_key, prop_label in _PROPERTIES_TO_COMPARE:
        row: dict[str, Any] = {"Property": prop_label}
        chart_data[prop_label] = []
        for name, compound in zip(names, compounds, strict=False):
            val = _get_property(compound, prop_key)
            row[name] = f"{val:.2f}" if isinstance(val, float) else (str(val) if val is not None else "—")
            chart_data[prop_label].append(val if isinstance(val, (int, float)) else None)
        rows.append(row)

    df = pd.DataFrame(rows)
    st.dataframe(df, use_container_width=True, hide_index=True)

    # Radar chart overlay
    _render_radar_overlay(chart_data, names)


def _render_radar_overlay(
    chart_data: dict[str, list[float | None]],
    names: list[str],
) -> None:
    """Render radar chart overlay for all selected compounds."""
    try:
        import plotly.graph_objects as go
    except ImportError:
        st.caption("Install plotly for radar chart visualization.")
        return

    # Filter to properties with at least one non-None value (allow negative values like XLogP)
    valid_props = {
        k: v for k, v in chart_data.items()
        if any(x is not None for x in v)
    }
    if len(valid_props) < 3:
        return

    categories = list(valid_props.keys())
    fig = go.Figure()

    for i, name in enumerate(names):
        values = [valid_props[prop][i] or 0 for prop in categories]
        values += values[:1]  # Close the polygon
        fig.add_trace(go.Scatterpolar(
            r=values,
            theta=categories + [categories[0]],
            name=name,
            fill="toself",
            opacity=0.4,
        ))

    fig.update_layout(
        polar={"radialaxis": {"visible": True}},
        title="Property Radar Chart",
        height=450,
    )
    st.plotly_chart(fig, use_container_width=True)


def _render_druglikeness_comparison(
    compounds: list[dict[str, Any]],
    names: list[str],
) -> None:
    """Render drug-likeness pass/fail comparison table."""
    st.subheader("Drug-likeness Comparison")

    rows: list[dict[str, Any]] = []
    for filter_name in _DRUGLIKENESS_FILTERS:
        row: dict[str, Any] = {"Filter": filter_name}
        for name, compound in zip(names, compounds, strict=False):
            passed = _get_druglikeness_result(compound, filter_name)
            if passed is None:
                row[name] = "—"
            elif passed:
                row[name] = "PASS"
            else:
                row[name] = "FAIL"
        rows.append(row)

    df = pd.DataFrame(rows)
    st.dataframe(df, use_container_width=True, hide_index=True)


def _render_admet_comparison(
    compounds: list[dict[str, Any]],
    names: list[str],
) -> None:
    """Render ADMET comparison if available for selected compounds."""
    # Check if any compound has ADMET data
    admet_keys: set[str] = set()
    for compound in compounds:
        dl = compound.get("druglikeness")
        if isinstance(dl, dict):
            admet_keys.update(k for k in dl if k.startswith("admet_"))

    if not admet_keys:
        with st.expander("ADMET Comparison"):
            st.caption("Run ADMET predictions on compounds first (via Batch Screen or Profile page).")
        return

    with st.expander("ADMET Comparison"):
        rows = []
        for key in sorted(admet_keys):
            row: dict[str, Any] = {"ADMET Property": key}
            for name, compound in zip(names, compounds, strict=False):
                dl = compound.get("druglikeness", {})
                val = dl.get(key) if isinstance(dl, dict) else None
                row[name] = f"{val:.3f}" if isinstance(val, float) else str(val) if val is not None else "—"
            rows.append(row)
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)


def _get_property(compound: dict[str, Any], key: str) -> Any:
    """Extract a property value from a compound dict, checking nested properties."""
    val = compound.get(key)
    if val is None:
        props = compound.get("properties", {})
        if isinstance(props, dict):
            val = props.get(key)
    if val is None:
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return val


def _get_druglikeness_result(compound: dict[str, Any], filter_name: str) -> bool | None:
    """Extract pass/fail result for a drug-likeness filter."""
    dl = compound.get("druglikeness")
    if dl is None:
        return None

    if isinstance(dl, dict):
        filters = dl.get("filters", {})
        if isinstance(filters, dict):
            result = filters.get(filter_name)
            if isinstance(result, dict):
                return result.get("pass")
            if isinstance(result, bool):
                return result

    if hasattr(dl, "filters"):
        filters = dl.filters
        if isinstance(filters, dict):
            result = filters.get(filter_name)
            if isinstance(result, dict):
                return result.get("pass")

    return None
